# Ported from iHEA_2027 R/functions/metrics.R at freeze-ihea2027-abstract-v3 (a79e102); changes: cvplus_interval renamed to .cvplus_interval, body unchanged; .conformal() (CV+ intervals under any value set, blocked by newdata rows) and CONFORMAL_BLOCK are new

# CV+ prediction interval (Barber et al. 2021). test_pred_by_fold: n_test x K
# matrix, column k the test predictions of the model fitted without fold k;
# oof_abs_resid: |R_i| of the n OOF rows; resid_fold: the fold (column index)
# of each OOF row. For each test person, lower = floor(alpha (n + 1))-th
# smallest of mu_{-k(i)}(x) - |R_i|, upper = ceiling((1 - alpha)(n + 1))-th
# smallest of mu_{-k(i)}(x) + |R_i|; -Inf / +Inf when the order is out of
# range (n too small). Not capped at 1 here. Returns list(lower, upper).
# A vector test_pred_by_fold is one test person (1 x K); the residuals are
# taken in absolute value whatever the caller passes.
.cvplus_interval <- function(test_pred_by_fold, oof_abs_resid, resid_fold, alpha) {
  tp <- if (is.null(dim(test_pred_by_fold))) matrix(test_pred_by_fold, nrow = 1) else as.matrix(test_pred_by_fold)
  oof_abs_resid <- abs(oof_abs_resid)
  n <- length(oof_abs_resid)
  if (length(resid_fold) != n) stop("cvplus_interval: resid_fold must have one fold per residual")
  if (!all(resid_fold %in% seq_len(ncol(tp)))) stop("cvplus_interval: resid_fold must index the columns of test_pred_by_fold")
  kl <- floor(alpha * (n + 1))
  ku <- ceiling((1 - alpha) * (n + 1))
  mu <- tp[, resid_fold, drop = FALSE]           # n_test x n
  r <- matrix(oof_abs_resid, nrow(tp), n, byrow = TRUE)
  lower <- if (kl < 1) rep(-Inf, nrow(tp)) else
    apply(mu - r, 1, function(v) sort(v, partial = kl)[kl])
  upper <- if (ku > n) rep(Inf, nrow(tp)) else
    apply(mu + r, 1, function(v) sort(v, partial = ku)[ku])
  list(lower = unname(lower), upper = unname(upper))
}

# Newdata rows per .cvplus_interval call (its work is n_test x n_train).
CONFORMAL_BLOCK <- 500L

# CV+ interval of level `level` for newdata under value set vs, from a bundle.
# Response: the point of fold k is mean_utility(fold k, newdata, vs); the
# out-of-fold prediction is the row mean of the OOF profiles scored under vs.
# Direct: the point of fold k is point(fold k, newdata) clamped to bundle$vs
# and, when a crosswalk is given, mapped to vs; the OOF prediction is oof$mu
# treated the same way. Residuals are |y - yhat| with y = .score(y_profiles,
# vs). Intervals are capped to [.vs_min(vs), 1].
.conformal <- function(bundle, newdata, vs, level, crosswalk = NULL) {
  if (!is.numeric(level) || length(level) != 1L || is.na(level) || level <= 0 || level >= 1) {
    stop("level must be a single number strictly between 0 and 1", call. = FALSE)
  }
  if (length(bundle$folds) == 0L) {
    stop("this bundle has no fold fits; refit with interval = TRUE", call. = FALSE)
  }
  m <- .p2e_method(bundle$method)
  direct <- identical(m$family, "direct")
  if (direct) {
    if (!inherits(bundle$vs, "p2e_vs")) stop("direct bundle has no value set", call. = FALSE)
    if (!inherits(vs, "p2e_vs")) {
      stop("a direct method needs the target vs to be a single value set", call. = FALSE)
    }
    if (!identical(vs$code, bundle$vs$code) && is.null(crosswalk)) {
      stop("the bundle was fitted on value set '", bundle$vs$code, "' but the target is '", vs$code,
           "'; a crosswalk is needed", call. = FALSE)
    }
  } else if (!is.null(crosswalk)) {
    stop("a crosswalk applies to direct methods only", call. = FALSE)
  }

  map_direct <- function(u) {
    u <- .clamp(u, bundle$vs)
    if (is.null(crosswalk)) u else .apply_crosswalk(u, crosswalk, vs)
  }

  y <- .score(bundle$y_profiles, vs)
  oof_hat <- if (direct) map_direct(bundle$oof$mu) else rowMeans(.score(bundle$oof$profiles, vs))
  resid <- abs(y - oof_hat)

  n <- nrow(newdata)
  tp <- matrix(vapply(bundle$folds, function(f) {
    if (direct) map_direct(m$point(f, newdata)) else m$mean_utility(f, newdata, vs)
  }, numeric(n)), n, length(bundle$folds))

  starts <- seq(1L, n, by = CONFORMAL_BLOCK)
  parts <- lapply(starts, function(s) {
    rows <- s:min(n, s + CONFORMAL_BLOCK - 1L)
    .cvplus_interval(tp[rows, , drop = FALSE], resid, bundle$fold_id, alpha = 1 - level)
  })
  lower <- unlist(lapply(parts, `[[`, "lower"))
  upper <- unlist(lapply(parts, `[[`, "upper"))
  list(lower = pmax(lower, .vs_min(vs)), upper = pmin(upper, 1))
}
