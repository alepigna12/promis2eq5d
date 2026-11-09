# Ported from iHEA_2027 R/functions/models/D12_klapproth.R at freeze-ihea2027-abstract-v3 (a79e102); changes: constructor .method_D12() instead of register_method(); formulas take baseenv() instead of globalenv(); stats functions namespaced; the draws are .d12_draw (keyed, with draw_type: as_published is mu + sigma z unclamped, corrected uses .censored_normal_location and clamps) instead of .d12_predict_draws; the fit stores no lm object

# D12: the Klapproth et al. (2020) specification, refitted (spec sec. 5.1, D12).
# The fair, leakage-free comparison with the published algorithm K-pub
# (klapproth2020.R): Klapproth's model search, rerun on each training set.
# One variant, list(form = "klapproth"), for every set:
#   candidates  OLS of y on the 23 predictors of .d12_formula(): AGE and
#               female, then T_x, I(T_x^2) and I(T_x^3) (raw, uncentred
#               powers) for each of the 7 T-scores. P9 is excluded, as in the
#               paper; the predictors are PREDICTOR_SETS$S3K (TSCORES + COVS).
#   selection   backward stepwise by BIC, step(direction = "backward",
#               k = log(n), trace = 0), with no hierarchy constraint (a
#               square or cube can stay without its lower powers). It runs
#               inside the fit, so it is repeated in every CV fold, bootstrap
#               replicate and final fit.
# The predictors are fixed by the specification: spec$set only labels the
# cell. The task grid creates S3K cells only (sets = "S3K").
# The fit stores terms_kept (term labels kept by step(), in formula order),
# n_candidates (23), beta and sigma; not the lm (it holds the training data).
# Point: X beta on the kept terms, unclamped (the cubes extrapolate beyond the
# training range). Draws: mu + sigma * z, with sigma the residual SD of the
# selected model. as_published keeps mu + sigma z unclamped (the engine
# clamps). corrected centres the draws on the clamped mean: the censored-normal
# location m with E[clamp(m + sigma Z)] = clamp(mu), then clamps the draws.
# Aliased coefficients (exact collinearity) get 0, as in D01.
# fit_info: AIC()/BIC() of the selected lm; n_par = its rank + 1 (sigma).

.method_D12 <- function() {
  list(
    id = "D12_klapproth", family = "direct",
    sets = METHOD_INFO$D12_klapproth$sets,
    variants = list(list(form = "klapproth")),
    aliases = c("D12", "klapproth"),
    needs = character(0),
    predict_needs = character(0),
    fit = function(train, set, variant, oof = NULL) {
      f <- .d12_fit(train, list(set = set, variant = variant))
      f$set <- set
      f$variant <- variant
      f
    },
    point = function(fit, newdata) .d12_predict_mean(fit, newdata),
    draw = function(fit, newdata, M, key, draw_type, lo) .d12_draw(fit, newdata, M, key, draw_type, lo)
  )
}

# The 23 candidate terms, in the order of the brief: COVS, then the powers of
# each T-score
.d12_candidates <- function() {
  c(COVS, as.vector(rbind(TSCORES, paste0("I(", TSCORES, "^2)"), paste0("I(", TSCORES, "^3)"))))
}

# Formula of the given terms (intercept-only without terms); baseenv(), so the
# formula keeps no reference to the caller's frame (the training data)
.d12_formula <- function(terms = .d12_candidates()) {
  rhs <- if (length(terms)) paste(terms, collapse = " + ") else "1"
  stats::as.formula(paste("y ~", rhs), env = baseenv())
}

.d12_fit <- function(train, spec) {
  df <- data.frame(y = train$y, train[PREDICTOR_SETS$S3K])
  if (anyNA(df)) stop("D12: training data has missing values")
  full <- stats::lm(.d12_formula(), data = df)
  # step() re-evaluates lm(formula, data = df) in this frame at every step
  sel <- stats::step(full, direction = "backward", k = log(nrow(df)), trace = 0)
  if (sel$df.residual < 1) stop("D12: no residual degrees of freedom")
  beta <- stats::coef(sel)
  beta[is.na(beta)] <- 0
  list(terms_kept = attr(stats::terms(sel), "term.labels"),
       n_candidates = length(attr(stats::terms(full), "term.labels")),
       beta = beta, sigma = summary(sel)$sigma,
       fit_info = list(converged = TRUE, aic = stats::AIC(sel), bic = stats::BIC(sel),
                       n_par = sel$rank + 1L))
}

.d12_predict_mean <- function(fit, newdata) {
  tt <- stats::delete.response(stats::terms(.d12_formula(fit$terms_kept)))
  X <- stats::model.matrix(tt, stats::model.frame(tt, newdata[PREDICTOR_SETS$S3K],
                                                  na.action = stats::na.pass))
  unname(drop(X[, names(fit$beta), drop = FALSE] %*% fit$beta))
}

# Keyed draws: K = max(key) units, row i uses unit key[i]. Rows sharing a key
# share their draws; key = seq_len(n) reproduces the pipeline's
# mu + matrix(rnorm(n * M, 0, sigma), n, M)
.d12_draw <- function(fit, newdata, M, key, draw_type, lo) {
  if (!identical(draw_type, "as_published") && !identical(draw_type, "corrected")) {
    stop("D12 draw_type must be \"as_published\" or \"corrected\", not ",
         paste(format(draw_type), collapse = " "), call. = FALSE)
  }
  mu <- .d12_predict_mean(fit, newdata)
  n <- length(mu)
  key <- as.integer(key)
  if (length(key) != n || anyNA(key) || min(key) < 1L) {
    stop("key must hold one positive integer per row of newdata", call. = FALSE)
  }
  K <- max(key)
  z <- matrix(stats::rnorm(K * M), K, M)[key, , drop = FALSE]
  sigma <- fit$sigma
  if (identical(draw_type, "as_published")) return(mu + sigma * z)
  target <- pmin(pmax(mu, lo), 1)
  m <- .censored_normal_location(target, sigma, lo, 1)
  pmin(pmax(m + sigma * z, lo), 1)
}
