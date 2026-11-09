# Ported from iHEA_2027 R/functions/models/R09_chain.R at freeze-ihea2027-abstract-v3 (a79e102); changes: constructor .method_R09() instead of register_method(); ML_THREADS replaced by the literal nthread = 1; ml_configs and the variant search not ported (variants[[1]] is a stand-in, Plan 2 supplies the selected variant); .r09_state_probs, .r09_profile_prob and R09_CHUNK not ported (R09_EXACT_MEAN = FALSE is in constants.R); fit() adds $set and $variant; each booster is stored as xgboost's portable raw model (ubj) and rebuilt for prediction (.r09_boosters), so a saved fit loads across xgboost versions and holds no external pointer; simulate() takes a key (K units, uniforms matrix(runif(K * M), K, M)[key, ]; key = seq_len(n) consumes the RNG as the pipeline); mean_utility() takes a p2e_vs value set and runs under withr::with_seed(MC_SEED) with pinned RNG kinds; marginals() added (level frequencies of the same N_MC_MEAN profiles); mean_utility() and marginals() process rows in blocks of R09_MEAN_BLOCK = 500 inside one MC_SEED stream (the Monte Carlo stream differs from an unblocked run for n > 500; simulate() is unchanged, and simulate with key = seq_len(n) is the pipeline-equivalent path); the conditional probabilities, chain order and sparse-level handling unchanged

# R09: classifier chain (xgboost, objective multi:softprob, num_class 5;
# spec sec. 5.2, R9). In the fixed order R09_ORDER, E4 -> E1 -> E3 -> E2 -> E5,
# the classifier of each dimension has the set's predictors plus the levels
# of the dimensions before it, entered as numeric predictors 1-5; in training
# they are the respondents' true levels.
#   P(state) = P(E4) P(E1 | E4) P(E3 | E4, E1) P(E2 | E4, E1, E3) P(E5 | ...).
# Predictors: the set's raw columns as a numeric matrix (.r09_x; female as
# 0/1). A missing, non-numeric or NA predictor is an error. fit$order is the
# chain order, fit$features the predictor names of each classifier (named by
# dimension). Every xgboost call uses nthread = 1; each booster's seed
# parameter is drawn from the R stream, in chain order.
# Sparse levels: labels 0-4 with num_class = 5; a level absent from train is
# set to 0 in every conditional and floored (floor_probs()). A dimension with
# a single observed level is an error.
# The fit holds, per dimension, the booster as a raw vector
# (xgboost::xgb.save.raw, "ubj") and the observed levels: no R booster, no
# environment, no training data. A booster loaded from raw has no thread
# setting (it would use every core), so .r09_boosters() sets nthread = 1.
# Rebuilding a 6 MB booster costs about 0.03 s, small against a prediction,
# so the five are rebuilt per call and not cached.
# Predictions: mean_utility and marginals use the Monte Carlo mean of
# N_MC_MEAN simulated profiles (blocks of R09_MEAN_BLOCK rows). simulate is sequential: in chain order, one
# n x M block of uniforms per dimension, each draw cut at the conditional of
# its own sampled earlier levels (predicted once per distinct person x
# earlier levels).

R09_ORDER <- c("E4", "E1", "E3", "E2", "E5")

.method_R09 <- function() {
  list(
    id = "R09_chain", family = "response",
    sets = METHOD_INFO$R09_chain$sets,
    variants = list(list(eta = 0.3, max_depth = 2, subsample = 1, colsample_bytree = 1,
                         min_child_weight = 1, nrounds = 10)),
    aliases = c("R09", "chain"),
    needs = "xgboost",
    predict_needs = "xgboost",
    fit = function(train, set, variant, oof = NULL) {
      f <- .r09_fit(train, list(set = set, variant = variant))
      f$set <- set
      f$variant <- variant
      f
    },
    simulate = function(fit, newdata, M, key) .r09_simulate(fit, newdata, M, key),
    mean_utility = function(fit, newdata, vs) .r09_mean_utility(fit, newdata, vs),
    marginals = function(fit, newdata) .r09_marginals(fit, newdata)
  )
}

# Rows per block of the Monte Carlo mean (mirrors MVOP_MEAN_BLOCK of R06):
# n x N_MC_MEAN x 5 at once would be gigabytes for thousands of rows
R09_MEAN_BLOCK <- 500L

# Monte Carlo mean utility: blocks of rows in one MC_SEED stream, keeping only
# the row means (key = seq_len(rows of the block))
.r09_mean_utility <- function(fit, newdata, vs, block = R09_MEAN_BLOCK) {
  if (!inherits(vs, "p2e_vs")) stop("vs must be a p2e_vs value set")
  n <- nrow(newdata)
  u <- numeric(n)
  .r09_mc(function() {
    for (b in split(seq_len(n), ceiling(seq_len(n) / block))) {
      prof <- .r09_simulate(fit, newdata[b, , drop = FALSE], N_MC_MEAN, seq_along(b))
      u[b] <<- rowMeans(matrix(.score(prof, vs), length(b)))
    }
  })
  u
}

# Per-dimension level frequencies of the same blocked Monte Carlo profiles
.r09_marginals <- function(fit, newdata, block = R09_MEAN_BLOCK) {
  n <- nrow(newdata)
  out <- lapply(DIMS, function(d) matrix(0, n, 5))
  .r09_mc(function() {
    for (b in split(seq_len(n), ceiling(seq_len(n) / block))) {
      prof <- .r09_simulate(fit, newdata[b, , drop = FALSE], N_MC_MEAN, seq_along(b))
      for (d in seq_along(DIMS)) {
        lv <- matrix(prof[, , d], length(b), N_MC_MEAN)
        for (l in 1:5) out[[d]][b, l] <<- rowMeans(lv == l)
      }
    }
  })
  stats::setNames(out, DIMS)
}

# Runs f under the fixed Monte Carlo seed with pinned RNG kinds, leaving the
# caller's stream alone
.r09_mc <- function(f) {
  .with_mc_seed(f())
}

# Raw predictors of a set as a numeric matrix with column names
.r09_x <- function(df, set) {
  cols <- PREDICTOR_SETS[[set]]
  miss <- setdiff(cols, names(df))
  if (length(miss)) stop("R09_chain: missing predictor(s): ", paste(miss, collapse = ", "))
  bad <- cols[!vapply(df[cols], is.numeric, NA)]                 # no silent as.numeric() of factor codes
  if (length(bad)) stop("R09_chain: predictor(s) not numeric: ", paste(bad, collapse = ", "))
  X <- matrix(vapply(cols, function(v) as.numeric(df[[v]]), numeric(nrow(df))), nrow(df), dimnames = list(NULL, cols))
  if (anyNA(X)) stop("R09_chain: missing values in the predictors")
  X
}

.r09_need_xgboost <- function() {
  if (!requireNamespace("xgboost", quietly = TRUE)) {
    stop("R09_chain needs the 'xgboost' package; install it with install.packages(\"xgboost\")")
  }
}

# The fit's boosters rebuilt from their raw models, one thread each (named by dimension)
.r09_boosters <- function(fit) {
  .r09_need_xgboost()
  lapply(fit$dims, function(d) {
    b <- xgboost::xgb.load.raw(d$raw)
    xgboost::`xgb.model.parameters<-`(b, list(nthread = 1))
    b
  })
}

# Conditional level probabilities of chain step k (n x 5) for the predictor
# rows X and the earlier levels prev (n x (k - 1), chain order): absent
# levels set to 0, floored. boosters: .r09_boosters(fit)
.r09_cond <- function(fit, k, X, prev, boosters = .r09_boosters(fit)) {
  d <- fit$order[k]
  Z <- cbind(X, matrix(as.numeric(prev), nrow(X), k - 1L, dimnames = list(NULL, fit$order[seq_len(k - 1L)])))
  q <- stats::predict(boosters[[d]], xgboost::xgb.DMatrix(Z[, fit$features[[d]], drop = FALSE], nthread = 1))
  if (is.null(dim(q))) q <- matrix(q, nrow(Z), 5, byrow = TRUE)   # row-major vector (older xgboost)
  q <- unname(q)
  q[, setdiff(1:5, fit$dims[[d]]$levels)] <- 0
  floor_probs(q)
}

.r09_fit <- function(train, spec) {
  .r09_need_xgboost()
  v <- spec$variant
  X <- .r09_x(train, spec$set)
  f <- list(set = spec$set, order = R09_ORDER, features = list(), dims = list())
  for (k in seq_along(R09_ORDER)) {
    d <- R09_ORDER[k]
    lev <- sort(unique(train[[d]]))
    if (length(lev) < 2) stop("R09_chain: only one level of ", d, " observed in train")
    f$features[[d]] <- c(colnames(X), R09_ORDER[seq_len(k - 1L)])
    Z <- cbind(X, as.matrix(train[R09_ORDER[seq_len(k - 1L)]]) + 0)   # true earlier levels, numeric
    params <- list(objective = "multi:softprob", num_class = 5L, eta = v$eta, max_depth = v$max_depth,
                   subsample = v$subsample, colsample_bytree = v$colsample_bytree,
                   min_child_weight = v$min_child_weight, nthread = 1,
                   seed = sample.int(.Machine$integer.max, 1))
    dtrain <- xgboost::xgb.DMatrix(Z[, f$features[[d]], drop = FALSE], label = train[[d]] - 1L, nthread = 1)
    booster <- xgboost::xgb.train(params = params, data = dtrain, nrounds = v$nrounds, verbose = 0)
    f$dims[[d]] <- list(raw = xgboost::xgb.save.raw(booster, raw_format = "ubj"), levels = lev)
  }
  bs <- .r09_boosters(f)
  obs <- as.matrix(train[R09_ORDER])
  ok <- vapply(seq_along(R09_ORDER), function(k) {
    all(is.finite(.r09_cond(f, k, X, obs[, seq_len(k - 1L), drop = FALSE], bs)))
  }, NA)
  f$fit_info <- list(converged = all(ok), aic = NA_real_, bic = NA_real_, n_par = NA_integer_)
  f
}

# Sequential simulation: integer n x M x 5 (DIMS order). Uniforms for K units
# (K = max(key)), one K x M block per chain step, row i using unit key[i]
.r09_simulate <- function(fit, newdata, M, key = seq_len(nrow(newdata))) {
  X <- .r09_x(newdata, fit$set)
  n <- nrow(X)
  key <- as.integer(key)
  if (length(key) != n || anyNA(key) || min(key) < 1L) stop("key must hold one positive integer per row of newdata")
  K <- max(key)
  bs <- .r09_boosters(fit)
  lv <- array(0L, c(n, M, 5))                       # chain order
  for (k in seq_along(fit$order)) {
    u <- matrix(stats::runif(K * M), K, M)[key, , drop = FALSE]
    prev <- matrix(lv[, , seq_len(k - 1L)], n * M, k - 1L)
    ckey <- (rep(seq_len(n), M) - 1) * 5^(k - 1) + drop((prev - 1) %*% 5^(seq_len(k - 1L) - 1))   # person x earlier levels
    first <- !duplicated(ckey)
    p <- .r09_cond(fit, k, X[rep(seq_len(n), M)[first], , drop = FALSE], prev[first, , drop = FALSE], bs)
    cp <- t(apply(p, 1, cumsum))[match(ckey, ckey[first]), , drop = FALSE]
    l <- rep(1L, n * M)
    for (j in 1:4) l <- l + (as.vector(u) > cp[, j])
    lv[, , k] <- l
  }
  out <- lv[, , match(DIMS, fit$order), drop = FALSE]
  storage.mode(out) <- "integer"
  out
}
