# Ported from iHEA_2027 R/functions/selection.R at freeze-ihea2027-abstract-v3 (a79e102); changes: .select_cell (lines 212-241) becomes .select_variant() on a plain table (variant_id, ok, mae, tie_metric), keeping SELECT_MAE_TIE and .TIE_EPS (lines 56, 60): candidates are the ok rows (the pipeline's usable units, lines 18-30), an NA MAE among them is an error, the tie metric is whatever the caller computed (-logp or CRPS) with no CRPS fallback, and only n_candidates and n_tied are added to the winning row; the cv stage's per-unit metrics (R/functions/tasks.R lines 90-100, 373-396, 667-690; selection.R .oof_metrics lines 94-98) become .cv_unit(); .crps_sample from R/functions/metrics.R lines 55-60 (crps_draws, used by tasks.R line 381 as scoringRules::crps_sample) written out in base R (the empirical-CDF formula; scoringRules is not a dependency); M_DRAWS_SELECT is tasks.R's M_DRAWS (line 134)

# Variant selection for the training run: .cv_unit() gives the out-of-fold
# metrics of one cell (method x set x target) x variant; .select_variant()
# picks a cell's variant from those units as the pipeline did.

SELECT_MAE_TIE <- 0.0005
.TIE_EPS <- sqrt(.Machine$double.eps)     # floating-point slack on the tie threshold (pipeline)
M_DRAWS_SELECT <- 200L                    # predictive draws per row for the CRPS (the pipeline's M_DRAWS)

# Per-row CRPS of predictive draws (n x M, row i = person i) against y: the
# CRPS of the empirical distribution of the M draws,
#   (2 / M) sum_j (1{y < x_(j)} - (j - 1/2) / M) (x_(j) - y),
# with x_(1) <= ... <= x_(M) the sorted draws of the row (the "edf" method of
# scoringRules::crps_sample, which the pipeline calls).
.crps_sample <- function(y, draws) {
  draws <- as.matrix(draws)
  n <- nrow(draws); M <- ncol(draws)
  if (length(y) != n) stop("y must have one value per row of draws")
  if (anyNA(y) || anyNA(draws)) stop("missing values in y or draws")
  X <- matrix(t(apply(draws, 1, sort)), n, M)
  a <- matrix((seq_len(M) - 0.5) / M, n, M, byrow = TRUE)
  2 / M * rowSums(((y < X) - a) * (X - y))
}

# The variant of one cell from its units: data.frame with variant_id, ok, mae
# and tie_metric, one row per variant. Candidates are the ok rows (a variant
# that failed in any fold is never chosen). Every candidate with mae within
# SELECT_MAE_TIE of the minimum (plus .TIE_EPS of floating-point slack) is
# tied; among the tied the lowest tie_metric wins (NA last), then the lowest
# mae, then the lowest variant_id. Without a tie, the minimum-mae candidate
# (lowest variant_id among equal mae). cell: an optional label (e.g.
# "D03_clad S3 UK") for the errors. Returns the winning row with n_candidates
# (ok rows) and n_tied.
.select_variant <- function(units, cell = NULL) {
  lab <- if (is.null(cell)) "" else paste0(" for ", cell)
  miss <- setdiff(c("variant_id", "ok", "mae", "tie_metric"), names(units))
  if (length(miss)) stop("units is missing column(s): ", paste(miss, collapse = ", "), call. = FALSE)
  d <- units[units$ok %in% TRUE, , drop = FALSE]
  if (!nrow(d)) stop("no usable variant", lab, ": every variant failed", call. = FALSE)
  if (anyNA(d$mae)) stop("a usable variant has a missing CV MAE", lab, call. = FALSE)
  best <- min(d$mae)
  if (!is.finite(best)) stop("non-finite CV MAE", lab, call. = FALSE)
  pick <- order(d$mae, d$variant_id)[1]                       # the minimum-MAE variant
  tied <- which(d$mae - best <= SELECT_MAE_TIE + .TIE_EPS)
  if (length(tied) > 1) {
    t <- d[tied, , drop = FALSE]
    pick <- tied[order(t$tie_metric, t$mae, t$variant_id, na.last = TRUE)[1]]
  }
  s <- d[pick, , drop = FALSE]
  s$n_candidates <- nrow(d)
  s$n_tied <- length(tied)
  rownames(s) <- NULL
  s
}

# Fit k of a cv unit: method$fit on the rows with fold_id != k under the
# bundle's fold seed (seed + 1000 + k), as .fit_bundle() fits fold k. train
# must already hold y (direct methods). fit() gets oof = NULL, so a uses_oof
# method (D03) also computes its own residual pairs over the remaining folds
# (the pipeline's cv fit; its core fit is the one the bundle makes with oof =
# FALSE, which runs first in the same stream). Returns .try_fit()'s
# list(fit, reason).
.cv_fold_fit <- function(method, train, set, variant, fold_id, k, seed) {
  train[["fold"]] <- fold_id
  .with_fit_seed(seed + 1000L + k, .try_fit(method, train[fold_id != k, , drop = FALSE], set, variant, NULL))
}

# Out-of-fold metrics of one variant of a cell (the pipeline's cv unit):
# list(ok, reason, mae, tie_metric, n). The K = max(fold_id) fold fits are
# .cv_fold_fit(); the first fold whose fit errors, does not converge or whose
# predictions error makes the unit not usable (ok = FALSE, reason "fold k:
# <message>", mae and tie_metric NA).
#   response  y = each row's utility on its own country's value set
#             (p2e_value_set(country)); mae = mean |y - mean_utility()| of
#             the fold fits on the same value sets; tie_metric = mean -logp
#             of the observed profile (NA for a method without logp)
#   direct    y = utility on p2e_value_set(target) for every row (the
#             rescored arm); mae = mean |y - point()| with the points clamped
#             to the value set; tie_metric = mean CRPS of M_DRAWS_SELECT
#             as_published draws per row, clamped (the pipeline's
#             predict_model draws), made under seed + 2000 + k (the bundle's
#             OOF seed); a uses_oof method draws from its fold fit's own pool
# n is the number of out-of-fold rows, nrow(train).
.cv_unit <- function(method, train, set, target, variant, fold_id, seed) {
  .require_pkgs(method$needs, paste0("fit ", method$id))
  n <- nrow(train)
  if (length(fold_id) != n || anyNA(fold_id) || any(!fold_id %in% seq_len(max(fold_id, 1L)))) {
    stop("fold_id must have length nrow(train) and values in 1..K", call. = FALSE)
  }
  K <- max(fold_id)
  if (any(tabulate(fold_id, K) == 0L)) stop("fold_id must use every fold 1..K (a fold is empty)", call. = FALSE)
  direct <- identical(method$family, "direct")
  prof <- as.matrix(train[DIMS])
  if (direct) {
    vs <- p2e_value_set(target)
    y <- .score(prof, vs)
    train[["y"]] <- y
  } else {
    cc <- as.character(train[["country"]])
    if (!length(cc) || anyNA(cc)) stop("train needs a country column without missing values", call. = FALSE)
    sets <- lapply(stats::setNames(unique(cc), unique(cc)), p2e_value_set)
    y <- .score(prof, unname(sets[cc]))
  }
  mu <- numeric(n)
  tie <- numeric(n)
  fail <- function(k, reason) {
    list(ok = FALSE, reason = paste0("fold ", k, ": ", reason), mae = NA_real_, tie_metric = NA_real_, n = n)
  }
  for (k in seq_len(K)) {
    r <- .cv_fold_fit(method, train, set, variant, fold_id, k, seed)
    if (is.null(r$fit)) return(fail(k, r$reason))
    te <- which(fold_id == k)
    nd <- train[te, , drop = FALSE]
    p <- tryCatch({
      if (direct) {
        m_k <- .clamp(method$point(r$fit, nd), vs)
        D <- .with_fit_seed(seed + 2000L + k,
                            method$draw(r$fit, nd, M_DRAWS_SELECT, seq_along(te), "as_published", .vs_min(vs)))
        list(mu = m_k, tie = .crps_sample(y[te], .clamp(D, vs)))
      } else {
        m_k <- numeric(length(te))
        for (g in unique(cc[te])) {
          w <- which(cc[te] == g)
          m_k[w] <- method$mean_utility(r$fit, nd[w, , drop = FALSE], sets[[g]])
        }
        list(mu = m_k, tie = if (is.null(method$logp)) rep(NA_real_, length(te)) else -method$logp(r$fit, nd))
      }
    }, error = function(e) conditionMessage(e))
    if (is.character(p)) return(fail(k, paste("prediction failed:", p)))
    mu[te] <- p$mu
    tie[te] <- p$tie
  }
  list(ok = TRUE, reason = NA_character_, mae = mean(abs(y - mu)), tie_metric = mean(tie), n = n)
}
