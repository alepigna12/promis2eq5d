# Ported from iHEA_2027 R/functions/design.R at freeze-ihea2027-abstract-v3 (a79e102); changes: dropped load_model_data, target_utility, state_weights, make_cv_folds, loco_folds, make_cohorts and PF_BAND_CUTS; lasso_select keeps the gaussian branch only (other families error) and checks that glmnet is installed; added stats:: and splines:: namespacing; pipeline sourcing and library() calls removed; PREDICTOR_SETS, ITEMS and TSCORES come from constants.R

# Design matrices and lasso selection: the data layer shared by every mapping model.
#
# Leakage rule: anything estimated from data (level merges, spline knots,
# centring, selected variables) is learned from the data given to
# build_design() when it is the training data, and is stored in `ref`. New data
# are only ever transformed with an existing `ref`.

MIN_LEVEL_N <- 20  # item levels with fewer respondents are merged (dummies form)


# Level merging for the "dummies" form. counts = respondents at levels 1..5.
# Rule: while some block of levels has fewer than `min_n` respondents, take the
# sparsest block (first one on ties) and merge it into the adjacent block
# towards the middle level (3): a block entirely below 3 merges upwards, a block
# entirely above 3 merges downwards. A block containing level 3 merges into its
# nearest neighbour, the one with fewer respondents (the lower one on ties).
# Stops when every block has at least `min_n` or only one block is left.
# Returns, for each level, the lowest level of its block.
.merge_levels <- function(counts, min_n = MIN_LEVEL_N) {
  blocks <- as.list(1:5)
  repeat {
    cnt <- vapply(blocks, function(b) sum(counts[b]), numeric(1))
    if (length(blocks) == 1 || all(cnt >= min_n)) break
    i <- which.min(cnt)
    k <- length(blocks)
    if (max(blocks[[i]]) < 3) {
      j <- i + 1L
    } else if (min(blocks[[i]]) > 3) {
      j <- i - 1L
    } else if (i == 1L) {
      j <- 2L
    } else if (i == k) {
      j <- i - 1L
    } else {
      j <- if (cnt[i + 1L] < cnt[i - 1L]) i + 1L else i - 1L
    }
    lo <- min(i, j)
    blocks[[lo]] <- sort(c(blocks[[i]], blocks[[j]]))
    blocks[[max(i, j)]] <- NULL
  }
  map <- integer(5)
  for (b in blocks) map[b] <- min(b)
  map
}

DESIGN_FORMS <- c("numeric", "dummies", "linear", "quadratic", "spline3", "spline3_int")

# Smooth terms: the T-scores and AGE present in the predictor set
.smooth_vars <- function(vars) intersect(vars, c(TSCORES, "AGE"))

.check_design_input <- function(df, vars) {
  miss <- setdiff(vars, names(df))
  if (length(miss)) stop("missing columns: ", paste(miss, collapse = ", "))
}

# Knots for a natural spline with df = 3: interior knots at the 1/3 and 2/3
# quantiles, boundary knots at the range. T-scores have ceilings/floors with
# many ties (for example a third of respondents at the maximum T_PF), so a
# quantile can sit on a boundary and ns() would fail or move it. Then the
# quantiles are taken among the values strictly inside the range, and
# duplicated knots are dropped; the spline then has fewer than 3 columns.
.spline_knots <- function(x) {
  bk <- range(x)
  k <- stats::quantile(x, c(1/3, 2/3), names = FALSE)
  if (!(all(k > bk[1] & k < bk[2]) && k[1] < k[2])) {
    inner <- x[x > bk[1] & x < bk[2]]
    k <- if (length(inner)) unique(stats::quantile(inner, c(1/3, 2/3), names = FALSE)) else numeric(0)
  }
  list(knots = k, boundary = bk)
}

# Learn everything from `df` (the reference/training data)
.learn_design <- function(df, set, form) {
  if (!set %in% names(PREDICTOR_SETS)) stop("unknown predictor set: ", set)
  if (!form %in% DESIGN_FORMS) stop("unknown form: ", form)
  if (form == "spline3_int" && set != "S3") stop("spline3_int is defined for S3 only")
  vars <- PREDICTOR_SETS[[set]]
  .check_design_input(df, vars)
  ref <- list(set = set, form = form, vars = vars)
  sv <- .smooth_vars(vars)
  if (form == "dummies") {
    items <- intersect(vars, ITEMS)
    ref$levmap <- lapply(stats::setNames(items, items), function(v) {
      x <- df[[v]]
      if (anyNA(x) || !all(x %in% 1:5)) stop("item ", v, " must be integers 1-5 without NA")
      .merge_levels(tabulate(x, 5))
    })
  }
  if (form %in% c("quadratic", "spline3_int")) {
    ref$center <- vapply(sv, function(v) mean(df[[v]]), numeric(1))
  }
  if (form %in% c("spline3", "spline3_int")) {
    ref$knots <- lapply(stats::setNames(sv, sv), function(v) .spline_knots(df[[v]]))
  }
  ref
}

# Apply `ref` to df: columns and values depend on `ref` only
.apply_design <- function(df, ref) {
  vars <- ref$vars; form <- ref$form
  .check_design_input(df, vars)
  cols <- list()
  sv <- .smooth_vars(vars)
  for (v in vars) {
    x <- df[[v]]
    if (anyNA(x)) stop("missing values in ", v)
    if (form == "dummies" && v %in% ITEMS) {
      if (!all(x %in% 1:5)) stop("item ", v, " must be integers 1-5")
      g <- ref$levmap[[v]][x]
      for (r in setdiff(sort(unique(ref$levmap[[v]])), min(ref$levmap[[v]]))) {
        cols[[paste0(v, "_", r)]] <- as.numeric(g == r)
      }
    } else if (form %in% c("spline3", "spline3_int") && v %in% sv) {
      k <- ref$knots[[v]]
      b <- splines::ns(x, knots = k$knots, Boundary.knots = k$boundary)
      for (j in seq_len(ncol(b))) cols[[paste0(v, "_ns", j)]] <- as.numeric(b[, j])
    } else {
      cols[[v]] <- as.numeric(x)
    }
  }
  if (form == "quadratic") {
    for (v in sv) cols[[paste0(v, "_sq")]] <- (df[[v]] - ref$center[[v]])^2
  }
  if (form == "spline3_int") {
    for (v in c("T_PF", "T_PI")) {
      cols[[paste0(v, "_x_AGE")]] <- (df[[v]] - ref$center[[v]]) * (df$AGE - ref$center[["AGE"]])
    }
  }
  X <- do.call(cbind, cols)
  rownames(X) <- NULL
  if (!all(is.finite(X))) stop("design matrix has non-finite values")
  X
}

# Design matrix for predictor set `set` ("S1", "S2", "S3", "S4", "S1D", "S3D")
# in one of the forms numeric, dummies, linear, quadratic, spline3, spline3_int.
#   numeric      predictors as they are
#   dummies      items as dummy variables; the lowest remaining level is the
#                reference; sparse levels (< 20 respondents in the ref data)
#                are merged, see .merge_levels()
#   linear       same columns as numeric (for T-score based sets)
#   quadratic    linear + squared (centred) T-scores and AGE
#   spline3      natural splines, df = 3, on each T-score and AGE; other
#                predictors enter linearly. Knots are stored in `ref`; new
#                values outside the boundary knots are extrapolated linearly
#   spline3_int  spline3 + T_PF x AGE and T_PI x AGE (centred); S3 only
# ref = NULL: everything is learned from `df`. Otherwise `ref` is applied
# unchanged to `df`. Returns list(X, ref).
build_design <- function(df, set, form, ref = NULL) {
  if (is.null(ref)) {
    ref <- .learn_design(df, set, form)
  } else if (ref$set != set || ref$form != form) {
    stop("ref was built for set ", ref$set, ", form ", ref$form)
  }
  list(X = .apply_design(df, ref), ref = ref)
}

# ---- lasso selection --------------------------------------------------------

# Names of the columns of X selected at lambda.1se (largest lambda whose CV
# error is within one standard error of the best), with CV folds `foldid`.
# Gaussian only: glmnet::cv.glmnet, lambda.1se.
lasso_select <- function(X, y, family = "gaussian", foldid) {
  if (!identical(family, "gaussian")) {
    stop("lasso_select supports family = \"gaussian\" only", call. = FALSE)
  }
  if (!requireNamespace("glmnet", quietly = TRUE)) {
    stop("Package 'glmnet' is needed for lasso_select(); install it with install.packages('glmnet')",
         call. = FALSE)
  }
  foldid <- as.integer(factor(foldid))  # character fold ids (e.g. country codes)
  X <- as.matrix(X)
  if (is.null(colnames(X))) stop("X needs column names")
  cv <- glmnet::cv.glmnet(X, y, family = family, foldid = foldid)
  b <- as.matrix(stats::coef(cv, s = "lambda.1se"))[-1, 1]
  names(b)[b != 0]
}
