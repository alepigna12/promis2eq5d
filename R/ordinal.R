# Ported from iHEA_2027 R/functions/model_interface.R (lines 480-680) at freeze-ihea2027-abstract-v3 (a79e102); changes: only floor_probs, .ordinal_design_form, .ordinal_conceptual_cols, .ordinal_location, .ordinal_clm and ordinal_marginals are ported; is_item_set added from the same file (R/functions/model_interface.R); .ordinal_location without the lasso_dummies branch (conceptual and the D01 forms only); stats functions imported via @importFrom

# Ordered-response helpers shared by the ordered-probit response models (R06).

#' @importFrom stats rnorm pnorm qnorm dnorm plogis dlogis optim optimize reformulate sd setNames
NULL

# Does the predictor set use the PROMIS items (S1, S2) or the T-scores?
is_item_set <- function(set) {
  if (!set %in% names(PREDICTOR_SETS)) stop("unknown predictor set: ", set)
  any(ITEMS %in% PREDICTOR_SETS[[set]])
}

# Floor every probability at eps, then renormalise each row to sum to 1
floor_probs <- function(p, eps = 1e-6) {
  p <- pmax(as.matrix(p), eps)
  p / rowSums(p)
}

# Conceptual variant: each dimension uses its CONCEPTUAL_SUBSETS predictors (item
# or T-score version), plus every predictor of the set that is not a PROMIS
# measure (COVS); items enter numerically, T-scores linearly.

.ordinal_design_form <- function(form, set) {
  switch(form,
         lasso_dummies = "dummies",
         conceptual    = if (is_item_set(set)) "numeric" else "linear",
         form)
}

.ordinal_conceptual_cols <- function(dim, set, colnames_x) {
  sub <- CONCEPTUAL_SUBSETS[[dim]][[if (is_item_set(set)) "items" else "tscores"]]
  miss <- setdiff(sub, colnames_x)
  if (length(miss)) stop("conceptual predictors not in the design: ", paste(miss, collapse = ", "))
  other <- setdiff(PREDICTOR_SETS[[set]], c(ITEMS, TSCORES, "P9"))
  intersect(colnames_x, c(sub, other))
}


# The location design of an ordered model: build_design() of the form's
# design form, and each dimension's columns (all or conceptual). label names the
# method in errors. Returns list(des, design_form, cols).
.ordinal_location <- function(train, spec, label) {
  form <- spec$variant$form
  if (!form %in% c(DESIGN_FORMS, "conceptual")) stop(label, ": unknown form ", form)
  design_form <- .ordinal_design_form(form, spec$set)
  des <- build_design(train, spec$set, design_form)
  cols <- setNames(lapply(DIMS, function(d) {
    if (form == "conceptual") .ordinal_conceptual_cols(d, spec$set, colnames(des$X))
    else colnames(des$X)
  }), DIMS)
  list(des = des, design_form = design_form, cols = cols)
}

# clm on standardised columns (means and SDs of this training data), then the
# coefficients are mapped back to the original scale: beta = beta_z / s and
# alpha = alpha_z + sum(beta_z * m / s). The model and its likelihood are
# unchanged; standardising only avoids clm's "nearly unidentifiable, rescale
# variables" failure (convergence code 3) when columns differ in scale by two
# orders of magnitude, as squared T-scores do. Converged = every clm code is 0.
# link: "logit" or "probit". scale_cols (R04): columns of X that also enter
# the scale, P(Y <= j) = F((alpha_j - x beta) / exp(z gamma)); they are
# standardised like the location, and mapped back as zeta = gamma_z / s and
# zeta0 = -sum(gamma_z * m / s), so the scale is exp(zeta0 + z zeta).
# Aliased coefficients (exact collinearity) get 0.
.ordinal_clm <- function(y, X, link, label, scale_cols = character(0)) {
  lev <- sort(unique(y))
  if (length(lev) < 2) stop(label, ": only one level observed in train")
  yf <- factor(y, levels = lev, ordered = TRUE)
  if (!ncol(X)) {
    fit <- ordinal::clm(y ~ 1, data = data.frame(y = yf), link = link)
    return(list(levels = lev, alpha = unname(fit$alpha), beta = numeric(0),
                converged = all(fit$convergence$code == 0), logLik = as.numeric(fit$logLik),
                edf = fit$edf, nobs = fit$nobs))
  }
  m <- colMeans(X)
  s <- apply(X, 2, sd)
  s[!is.finite(s) | s == 0] <- 1                       # constant column: aliased with the thresholds
  Z <- sweep(sweep(X, 2, m), 2, s, "/")
  dat <- data.frame(y = yf, Z, check.names = FALSE)
  fit <- if (!length(scale_cols)) ordinal::clm(y ~ ., data = dat, link = link) else {
    if (!all(scale_cols %in% colnames(X))) stop(label, ": scale columns not in the location design")
    ordinal::clm(y ~ ., scale = reformulate(paste0("`", scale_cols, "`")), data = dat, link = link)
  }
  bz <- unname(fit$beta)
  bz[is.na(bz)] <- 0
  out <- list(levels = lev, alpha = unname(fit$alpha) + sum(bz * m / s), beta = bz / s,
              converged = all(fit$convergence$code == 0), logLik = as.numeric(fit$logLik),
              edf = fit$edf, nobs = fit$nobs)
  if (length(scale_cols)) {
    gz <- unname(fit$zeta)
    gz[is.na(gz)] <- 0
    ms <- m[scale_cols]; ss <- s[scale_cols]
    out$scale_cols <- scale_cols
    out$zeta <- gz / ss
    out$zeta0 <- -sum(gz * ms / ss)
    out$converged <- out$converged && all(is.finite(c(out$alpha, out$beta, out$zeta)))
  }
  out
}


# Marginal probabilities of an ordered model: list over DIMS of n x 5
# matrices (levels 1-5), floored. P(Y <= j) = F(alpha_j - x beta), as in clm,
# divided by the scale exp(zeta0 + z zeta) for a dimension with one (R04).
ordinal_marginals <- function(fit, newdata) {
  X <- build_design(newdata, fit$ref$set, fit$design_form, fit$ref)$X
  n <- nrow(X)
  cdf <- switch(fit$link, logit = plogis, probit = pnorm, stop("unknown link: ", fit$link))
  setNames(lapply(DIMS, function(d) {
    m <- fit$dims[[d]]
    eta <- if (length(m$beta)) drop(X[, fit$cols[[d]], drop = FALSE] %*% m$beta) else rep(0, n)
    lin <- outer(-eta, m$alpha, "+")
    if (length(m$scale_cols)) lin <- lin / exp(m$zeta0 + drop(X[, m$scale_cols, drop = FALSE] %*% m$zeta))
    cp <- cdf(lin)
    p <- matrix(0, n, 5)
    p[, m$levels] <- cbind(cp, 1) - cbind(0, cp)
    floor_probs(p)
  }), DIMS)
}
