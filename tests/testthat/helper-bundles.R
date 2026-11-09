# Test helpers for fit bundles.

# A direct method for tests (not registered): OLS of train$y on the S4 columns.
# fit() counts its own calls (point, then K folds, then bootstrap b = call
# 1 + K + b) and errors on bootstrap replicate `fail_on`. With uses_oof = TRUE
# the point fit is made after the folds and the oof it receives is kept in
# fit$oof_given.
toy_method <- function(fail_on = NULL, uses_oof = FALSE, K = 3, bad_call = NULL, err_call = NULL, fi_bad_call = NULL) {
  calls <- new.env()
  calls$n <- 0L
  cols <- PREDICTOR_SETS$S4
  list(
    id = "toy", family = "direct", sets = "S4", variants = list(list()),
    aliases = "toy", needs = character(0), predict_needs = character(0),
    uses_oof = uses_oof,
    fit = function(train, set, variant, oof = NULL) {
      calls$n <- calls$n + 1L
      if (!is.null(err_call) && calls$n == err_call) stop("toy fit errored on purpose")
      if (!is.null(fail_on) && !uses_oof && calls$n == 1L + K + fail_on) {
        stop("toy fit failed on purpose")
      }
      # lm.fit on a plain matrix: the fit holds no environments, so it saves silently
      X <- cbind(1, as.matrix(train[cols]))
      m <- stats::lm.fit(X, train$y)
      out <- list(coef = unname(m$coefficients), sd = sqrt(sum(m$residuals^2) / max(1, nrow(X) - ncol(X))),
           set = set, variant = variant, oof_given = oof, ids = train$ID,
           n_train = nrow(train), draw_probe = stats::runif(1))
      if (!is.null(bad_call)) out$converged <- calls$n != bad_call
      if (!is.null(fi_bad_call)) out$fit_info <- list(converged = calls$n != fi_bad_call)
      out
    },
    point = function(fit, newdata) as.numeric(cbind(1, as.matrix(newdata[cols])) %*% fit$coef),
    draw = function(fit, newdata, M, key, draw_type, lo) {
      mu <- as.numeric(cbind(1, as.matrix(newdata[cols])) %*% fit$coef)
      mu + matrix(stats::rnorm(length(mu) * M, 0, fit$sd), length(mu), M)
    }
  )
}

.stand_in_env <- new.env()

# A cached bundle fitted on synthetic data (K = 3, B = 4, seed 1, first
# variant), written to a session temp directory that becomes promis2eq5d.fit_dir.
stand_in_bundle <- function(method, set, target = NA) {
  key <- paste(method, set, if (is.na(target)) "any" else target, sep = "_")
  if (is.null(.stand_in_env$dir)) {
    .stand_in_env$dir <- tempfile("p2e-fits-")
    dir.create(.stand_in_env$dir)
  }
  options(promis2eq5d.fit_dir = .stand_in_env$dir)
  if (is.null(.stand_in_env[[key]])) {
    m <- .p2e_method(method)
    train <- make_synthetic_data(400, seed = 1, countries = c("DE", "UK"))
    b <- .fit_bundle(m, train, set, target, .p2e_method_variants(m, set)[[1]], K = 3, B = 4, seed = 1)
    saveRDS(b, file.path(.stand_in_env$dir, paste0(key, ".rds")))
    .stand_in_env[[key]] <- b
  }
  .stand_in_env[[key]]
}

# Makes sure the stand-in bundles that predict_eq5d() tests need exist.
use_stand_in_fits <- function() {
  for (s in c("S1", "S2", "S3")) stand_in_bundle("R06_mvop_dim", s)
  stand_in_bundle("R09_chain", "S3")
  for (s in c("S3", "S1")) for (tg in c("FR", "DE", "UK")) stand_in_bundle("D03_clad", s, tg)
  for (tg in c("FR", "DE", "UK")) stand_in_bundle("D12_klapproth", "S3K", tg)
  invisible(TRUE)
}
