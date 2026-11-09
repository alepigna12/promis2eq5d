# Ported from iHEA_2027 R/functions/mvop.R at freeze-ihea2027-abstract-v3 (a79e102); changes: mvop_profile_prob, MVOP_ABSEPS, MVOP_SEED and register_mvop_response not ported (profile_prob is evaluation-only, so mvtnorm is not needed; pbivnorm and Matrix are called with ::); mvop_simulate_profiles takes a key (K units, row (i, m) uses unit key[i] + K * (m - 1); key = seq_len(n) consumes the RNG as the pipeline); mvop_mean_utility takes a p2e_vs value set via .vs_additive() and .clamp(), runs its Monte Carlo part under withr::with_seed(MC_SEED) with pinned RNG kinds; other MVOP_* constants and the estimator unchanged
# Shared code of the correlated ordered response models (spec §5.2, "Response
# models R5, R6 and R10: estimation"; plan Task 15): the multivariate ordered
# probit R05 (D01 forms) and R06 (conceptual predictors per dimension), and
# the Gaussian copula with ordered-logit margins R10 (R01's variants). In
# TASK_CODE_FILES, so it is hashed for every method.
#
#   f  <- mvop_fit(train, spec, link = "probit", label = "R05")   # or link = "logit" (R10)
#   mg <- ordinal_marginals(f, newdata)                          # the margins, floored
#   s  <- mvop_simulate_profiles(f, newdata, M)                  # n x M x 5
#   p  <- mvop_profile_prob(f, newdata, profiles)                # n
#   mu <- mvop_mean_utility(f, newdata, vt)                      # n
#
# Model. Each dimension j has a latent z_j and P(Y_j <= l | x) = F(alpha_jl -
# x_j' beta_j), F = pnorm (R05, R06) or plogis (R10), with x_j the
# dimension's columns of .ordinal_location() (all, its own lasso-selected or
# its conceptual columns). z is N(0, R) with a 5 x 5 correlation matrix R,
# constant across people. On the normal scale the cut points are c_ijl =
# alpha_jl - x_ij' beta_j (probit) or qnorm(plogis(alpha_jl - x_ij' beta_j))
# (logit: a Gaussian copula), clamped to +-MVOP_CLAMP; an infinite end is
# +-MVOP_BIG in the likelihood (pnorm(-10) is 8e-24).
#
# Objective: the negative pairwise log-likelihood (mvord's), the sum over the
# 10 pairs j < k and the people i of -log P(Y_ij, Y_ik), each term the
# bivariate normal probability (pbivnorm) of the rectangle (c_ij,y-1, c_ij,y]
# x (c_ik,y-1, c_ik,y] at rho_jk.
#
# Estimation: iterated two-stage (block coordinate descent), no multistart and
# no warm start.
#   Round 0: each margin alone with .ordinal_clm() (link probit: exactly R02;
#     logit: exactly R01), then each rho_jk by optimize() on (-MVOP_RHO_MAX,
#     MVOP_RHO_MAX) with the margins fixed: the two-stage estimator.
#   Each later round: (a) each dimension in turn, its thresholds and
#     coefficients by BFGS with the analytic gradient (correlations and the
#     other margins fixed), on the standardised design of .ordinal_clm() and
#     with ordered thresholds (alpha_1, log of the increments); (b) each rho_jk
#     as in round 0. A step that would worsen the objective is not taken, so
#     the objective never increases within a round (only a nearPD projection
#     can move it up).
#   Stop when a round gains less than MVOP_TOL, or after MVOP_MAX_ITER
#     rounds; reaching MVOP_MAX_ITER without meeting the tolerance keeps the
#     fit and sets fit_info$flag_max_iter (not a failure). MVOP_MAX_ITER = 0
#     is round 0 only, the plain two-stage estimator: the emergency switch for
#     a short cluster budget. It is a constant of this hashed file, so results
#     of the two settings are never mixed.
#   After each round, an R with smallest eigenvalue below MVOP_PD_EIG is
#     projected with Matrix::nearPD(corr = TRUE), counted in n_pd_fix.
# A coefficient that clm aliased (0) stays 0. The margins of round 0 are kept
# as clm returned them, so with MVOP_MAX_ITER = 0 they are R02's (R01's)
# numbers exactly; a margin refitted in a later round is mapped back to the
# original scale as in .ordinal_clm().
#
# Levels missing from train (as R01/R02): a margin has thresholds between
# its observed levels only, a missing level is never simulated, and it gets
# the floored probability in the margins (floor_probs()) and in profile_prob.
#
# Predictions:
#   simulate_profiles  z ~ N(0, R) via chol(R), cut at each person's c_ijl; the
#                      normals are drawn for K units (K = max(key)) and row
#                      (i, m) uses unit key[i] + K * (m - 1)
#   mean_utility       the exact additive part from the floored margins,
#                      intercept - sum_j sum_l P_ij(l) d(j, l), with the
#                      decrements of .vs_additive(), plus the mean over
#                      N_MC_MEAN simulated profiles of the value set's
#                      remainder (table - additive value); the simulation is
#                      skipped for a value set without remainder (DE, UK), and
#                      runs under the fixed seed MC_SEED otherwise
#   profile_prob       (pmvnorm over the profile's 5-d rectangle) is not ported:
#                      it is used for evaluation only, so mvtnorm is not needed
#
# fit object: as ordinal_fit() (ref, design_form, cols, dims, link; so
# ordinal_marginals() gives the margins), plus R, coef (a predictors x DIMS
# matrix on the original scale, 0 for a column a dimension does not use),
# objective_trace (the objective after rounds 0, 1, ...) and fit_info:
# converged, R, n_iter (rounds after round 0), flag_max_iter, n_pd_fix,
# bfgs_codes (optim()'s convergence code of every margin step, rounds x
# DIMS order; recorded only: a step is kept only if it lowers the
# objective, so a non-zero code does not make the fit fail), flag_bfgs_code
# (any non-zero code; counted in run_report.csv; author, 2026-11-02),
# objective, objective_round0 (negative pairwise log-likelihood), n_par
# (thresholds, coefficients not aliased, 10 correlations); aic and bic NA
# (a pairwise likelihood has no AIC).

MVOP_TOL      <- 1e-3    # stop when a round gains less than this in the objective
MVOP_MAX_ITER <- 10L     # rounds after round 0; 0 = plain two-stage (emergency switch)
MVOP_RHO_MAX  <- 0.99    # rho_jk searched on (-MVOP_RHO_MAX, MVOP_RHO_MAX)
MVOP_RHO_TOL  <- 1e-7    # optimize() tolerance for rho_jk
MVOP_PD_EIG   <- 1e-8    # smallest eigenvalue of R below this: nearPD projection
MVOP_CLAMP    <- 8       # normal-scale cut points of logit margins clamped to +-8
MVOP_BIG      <- 10      # an infinite end (and a probit cut point beyond it) in the likelihood

# ---- cut points ----------------------------------------------------------------

# Normal-scale cut points of linear cut values lin (any shape): probit lin
# itself, logit qnorm(plogis(lin)) computed on the log scale in the lower
# tail (so it stays finite), clamped to +-MVOP_CLAMP; then infinite or larger
# values at +-MVOP_BIG
.mvop_normal_cut <- function(lin, link) {
  out <- lin
  if (link == "logit") {
    fin <- is.finite(lin)
    a <- -abs(lin[fin])
    out[fin] <- -sign(lin[fin]) * qnorm(plogis(a, log.p = TRUE), log.p = TRUE)
    out[fin] <- pmin(pmax(out[fin], -MVOP_CLAMP), MVOP_CLAMP)
  } else if (link != "probit") stop("unknown link: ", link)
  out[] <- pmin(pmax(out, -MVOP_BIG), MVOP_BIG)
  out
}

# Derivative of the probability of a rectangle side with respect to its
# linear cut value: d Phi(c(lin)) / d lin = dnorm (probit) or dlogis (logit),
# 0 where the cut point is infinite or clamped
.mvop_cut_density <- function(lin, c, link) {
  d <- if (link == "logit") dlogis(lin) else dnorm(lin)
  lim <- if (link == "logit") MVOP_CLAMP else MVOP_BIG
  d[!is.finite(lin) | abs(c) >= lim] <- 0
  d
}

# Thresholds from their unconstrained parameters (alpha_1, log increments) and back
.mvop_g2a <- function(g) cumsum(c(g[1], exp(g[-1])))
.mvop_a2g <- function(a) c(a[1], log(diff(a)))

# One margin's rectangle for each person: theta = c(g, b) on the standardised
# design Z (n x length(b)), y the level index 1..K among the margin's
# observed levels. Returns list(lo, hi) of normal-scale cut points and
# list(dlo, dhi), the derivatives of Phi(lo), Phi(hi) by their linear values.
.mvop_rect <- function(theta, y, Z, link) {
  nt <- length(theta) - ncol(Z)
  A <- c(-Inf, .mvop_g2a(theta[seq_len(nt)]), Inf)
  eta <- if (ncol(Z)) drop(Z %*% theta[-seq_len(nt)]) else rep(0, length(y))
  lhi <- A[y + 1] - eta; llo <- A[y] - eta
  hi <- .mvop_normal_cut(lhi, link); lo <- .mvop_normal_cut(llo, link)
  list(lo = lo, hi = hi, dlo = .mvop_cut_density(llo, lo, link), dhi = .mvop_cut_density(lhi, hi, link))
}

# ---- pairwise likelihood -------------------------------------------------------

# Bivariate normal CDF at correlation r. A side at -MVOP_BIG gives 0 and a
# side at +MVOP_BIG the other side's pnorm (exact to 1e-23), so pbivnorm runs
# only on the finite pairs: most people are at level 1 of most dimensions.
.mvop_F2 <- function(x, y, r) {
  out <- numeric(length(x))
  xl <- x <= -MVOP_BIG; yl <- y <= -MVOP_BIG; xh <- x >= MVOP_BIG; yh <- y >= MVOP_BIG
  i <- xh & !yl; out[i] <- pnorm(y[i])
  i <- yh & !xl & !xh; out[i] <- pnorm(x[i])
  i <- !(xl | yl | xh | yh)
  if (any(i)) out[i] <- pbivnorm::pbivnorm(x[i], y[i], rep(r, sum(i)))
  out
}

# Probability of each person's pair of rectangles at correlation r
.mvop_pair_prob <- function(cj, ck, r) {
  p <- .mvop_F2(cj$hi, ck$hi, r) - .mvop_F2(cj$lo, ck$hi, r) - .mvop_F2(cj$hi, ck$lo, r) + .mvop_F2(cj$lo, ck$lo, r)
  pmax(p, 1e-300)
}

.mvop_pair_nll <- function(cj, ck, r) -sum(log(.mvop_pair_prob(cj, ck, r)))

# The objective: rects is a list of 5 .mvop_rect() results, R the correlations
.mvop_objective <- function(rects, R) {
  s <- 0
  for (j in 1:4) for (k in (j + 1):5) s <- s + .mvop_pair_nll(rects[[j]], rects[[k]], R[j, k])
  s
}

# The objective's terms that involve margin j (its 4 pairs), as a function of
# its parameters theta; m holds the margin's y, Z and link
.mvop_margin_nll <- function(theta, j, m, rects, R) {
  cj <- .mvop_rect(theta, m$y, m$Z, m$link)
  s <- 0
  for (k in setdiff(1:5, j)) s <- s + .mvop_pair_nll(cj, rects[[k]], R[j, k])
  s
}

# Analytic gradient of .mvop_margin_nll() in theta. For a pair, d P / d hi_j
# = dPhi(hi_j) * [Phi((hi_k - r hi_j) / s) - Phi((lo_k - r hi_j) / s)], s =
# sqrt(1 - r^2), and likewise (negated) for lo_j; hi_j moves with alpha_{j,y}
# and -x' b, lo_j with alpha_{j,y-1} and -x' b; then the chain rule to the
# threshold parameters (alpha_m = g_1 + sum_{i <= m} exp(g_i)).
.mvop_margin_grad <- function(theta, j, m, rects, R) {
  nt <- length(theta) - ncol(m$Z)
  g <- theta[seq_len(nt)]
  cj <- .mvop_rect(theta, m$y, m$Z, m$link)
  w <- numeric(length(m$y)); ga <- numeric(nt)
  for (k in setdiff(1:5, j)) {
    r <- R[j, k]; s <- sqrt(1 - r^2); ck <- rects[[k]]
    pr <- .mvop_pair_prob(cj, ck, r)
    dh <-  cj$dhi * (pnorm((ck$hi - r * cj$hi) / s) - pnorm((ck$lo - r * cj$hi) / s)) / pr
    dl <- -cj$dlo * (pnorm((ck$hi - r * cj$lo) / s) - pnorm((ck$lo - r * cj$lo) / s)) / pr
    w <- w + dh + dl
    for (l in seq_len(nt)) ga[l] <- ga[l] - sum(dh[m$y == l]) - sum(dl[m$y == l + 1])
  }
  gg <- vapply(seq_len(nt), function(i) if (i == 1) sum(ga) else exp(g[i]) * sum(ga[i:nt]), 0)
  c(gg, if (ncol(m$Z)) drop(crossprod(m$Z, w)))
}

# Each rho_jk by optimize() given the margins. R_old (later rounds): a rho
# that optimize() finds no better than the current one is kept.
.mvop_rho_step <- function(rects, R_old = NULL) {
  R <- diag(5)
  for (j in 1:4) for (k in (j + 1):5) {
    f <- function(r) .mvop_pair_nll(rects[[j]], rects[[k]], r)
    o <- optimize(f, c(-MVOP_RHO_MAX, MVOP_RHO_MAX), tol = MVOP_RHO_TOL)
    r <- o$minimum
    if (!is.null(R_old) && f(R_old[j, k]) <= o$objective) r <- R_old[j, k]
    R[j, k] <- R[k, j] <- r
  }
  R
}

# R if positive definite (smallest eigenvalue >= MVOP_PD_EIG), else its
# nearest correlation matrix; list(R, fixed)
.mvop_pd <- function(R) {
  if (min(eigen(R, symmetric = TRUE, only.values = TRUE)$values) >= MVOP_PD_EIG) return(list(R = R, fixed = FALSE))
  P <- as.matrix(Matrix::nearPD(R, corr = TRUE)$mat)
  dimnames(P) <- dimnames(R)
  list(R = P, fixed = TRUE)
}

# ---- fit -----------------------------------------------------------------------

# Fit R05/R06 (link "probit") or R10 (link "logit"); label names the method
# in errors. See the header for the estimator.
mvop_fit <- function(train, spec, link, label) {
  loc <- .ordinal_location(train, spec, label)
  X <- loc$des$X
  dims <- setNames(lapply(DIMS, function(d) {
    .ordinal_clm(train[[d]], X[, loc$cols[[d]], drop = FALSE], link, label)
  }), DIMS)
  # each margin on its standardised free columns (not aliased by clm), with
  # the round-0 parameters on that scale
  mg <- lapply(DIMS, function(d) {
    m0 <- dims[[d]]
    Xd <- X[, loc$cols[[d]], drop = FALSE]
    free <- which(m0$beta != 0)
    mu <- colMeans(Xd)[free]
    s <- apply(Xd, 2, sd)[free]
    s[!is.finite(s) | s == 0] <- 1
    Z <- sweep(sweep(Xd[, free, drop = FALSE], 2, mu), 2, s, "/")
    b <- m0$beta[free]
    list(y = match(train[[d]], m0$levels), Z = Z, link = link, free = free, mu = mu, s = s,
         theta = c(.mvop_a2g(m0$alpha - sum(b * mu)), b * s))
  })
  rects <- lapply(mg, function(m) .mvop_rect(m$theta, m$y, m$Z, link))
  pd <- .mvop_pd(.mvop_rho_step(rects))
  R <- pd$R; n_pd_fix <- as.integer(pd$fixed)
  obj <- .mvop_objective(rects, R)
  trace <- obj
  n_iter <- 0L; gain <- Inf; bfgs_codes <- integer(0)
  while (n_iter < MVOP_MAX_ITER) {
    for (j in 1:5) {
      m <- mg[[j]]
      o <- optim(m$theta, .mvop_margin_nll, .mvop_margin_grad, j = j, m = m, rects = rects, R = R,
                 method = "BFGS", control = list(maxit = 500, reltol = 1e-12))
      bfgs_codes <- c(bfgs_codes, as.integer(o$convergence))
      if (all(is.finite(o$par)) && is.finite(o$value) &&
          o$value < .mvop_margin_nll(m$theta, j, m, rects, R)) {
        mg[[j]]$theta <- o$par
        mg[[j]]$refit <- TRUE
        rects[[j]] <- .mvop_rect(o$par, m$y, m$Z, link)
      }
    }
    pd <- .mvop_pd(.mvop_rho_step(rects, R))
    R <- pd$R; n_pd_fix <- n_pd_fix + as.integer(pd$fixed)
    new <- .mvop_objective(rects, R)
    n_iter <- n_iter + 1L
    gain <- obj - new; obj <- new; trace <- c(trace, obj)
    if (gain < MVOP_TOL) break
  }
  # refitted margins back to the original scale, as in .ordinal_clm()
  for (j in which(vapply(mg, function(m) isTRUE(m$refit), NA))) {
    m <- mg[[j]]; nt <- length(m$theta) - length(m$free)
    bz <- m$theta[-seq_len(nt)]
    beta <- numeric(length(dims[[j]]$beta)); beta[m$free] <- bz / m$s
    dims[[j]]$alpha <- .mvop_g2a(m$theta[seq_len(nt)]) + sum(bz * m$mu / m$s)
    dims[[j]]$beta <- beta
  }
  dimnames(R) <- list(DIMS, DIMS)
  coef <- matrix(0, ncol(X), 5, dimnames = list(colnames(X), DIMS))
  for (d in DIMS) coef[loc$cols[[d]], d] <- dims[[d]]$beta
  n_par <- sum(vapply(mg, function(m) length(m$theta), 0)) + 10
  converged <- all(vapply(dims, `[[`, TRUE, "converged")) &&
    all(is.finite(unlist(lapply(dims, `[`, c("alpha", "beta"))))) && all(is.finite(R)) && is.finite(obj)
  list(ref = loc$des$ref, design_form = loc$design_form, cols = loc$cols, dims = dims, link = link,
       R = R, coef = coef, objective_trace = trace,
       fit_info = list(converged = converged, aic = NA_real_, bic = NA_real_, n_par = n_par, R = R,
                       n_iter = n_iter, flag_max_iter = MVOP_MAX_ITER > 0 && n_iter == MVOP_MAX_ITER && gain >= MVOP_TOL,
                       n_pd_fix = n_pd_fix, objective = obj, objective_round0 = trace[1],
                       bfgs_codes = bfgs_codes, flag_bfgs_code = any(bfgs_codes != 0L)))
}

# ---- predictions ---------------------------------------------------------------


# Normal-scale cut points of newdata: list over DIMS of n x (K_d - 1)
# matrices (the dimension's observed levels K_d)
.mvop_cuts <- function(fit, newdata) {
  X <- build_design(newdata, fit$ref$set, fit$design_form, fit$ref)$X
  n <- nrow(X)
  setNames(lapply(DIMS, function(d) {
    m <- fit$dims[[d]]
    eta <- if (length(m$beta)) drop(X[, fit$cols[[d]], drop = FALSE] %*% m$beta) else rep(0, n)
    .mvop_normal_cut(outer(-eta, m$alpha, "+"), fit$link)
  }), DIMS)
}

# Profiles: z ~ N(0, R) (K * M rows of rnorm through chol(R), K = max(key)); the
# level of a dimension is its observed level number 1 + #{l: z > c_l}. Row i,
# draw m uses unit row key[i] + K * (m - 1), so rows sharing a key share their
# normals; key = seq_len(n) is the pipeline's layout. Integer n x M x 5.
mvop_simulate_profiles <- function(fit, newdata, M, key = seq_len(nrow(newdata))) {
  cuts <- .mvop_cuts(fit, newdata)
  n <- nrow(newdata)
  key <- as.integer(key)
  if (length(key) != n || anyNA(key) || min(key) < 1L) stop("key must hold one positive integer per row of newdata")
  K <- max(key)
  Z <- matrix(rnorm(K * M * 5), K * M, 5) %*% chol(fit$R)
  unit <- rep(key, M) + K * rep(seq_len(M) - 1L, each = n)
  Z <- Z[unit, , drop = FALSE]
  out <- array(0L, c(n, M, 5))
  person <- rep(seq_len(n), M)
  for (d in seq_along(DIMS)) {
    cm <- cuts[[d]]
    k <- rep(1L, n * M)
    for (l in seq_len(ncol(cm))) k <- k + (Z[, d] > cm[person, l])
    out[, , d] <- fit$dims[[d]]$levels[k]
  }
  storage.mode(out) <- "integer"
  out
}

# Expected utility on the value set vs (one p2e_vs): the exact additive part
# from the floored margins plus the Monte Carlo mean (N_MC_MEAN profiles) of the
# remainder, skipped for value sets without one. The profiles are simulated
# in blocks of MVOP_MEAN_BLOCK rows (memory: n x N_MC_MEAN x 5 at once would
# be 400 MB for 5,000 rows) under withr::with_seed(MC_SEED), so the mean is
# deterministic and the caller's RNG stream is left alone. Clamped to
# [.vs_min(vs), 1] (the Monte Carlo part could cross it by about 1e-5 at most).
MVOP_MEAN_BLOCK <- 500L
mvop_mean_utility <- function(fit, newdata, vs) {
  if (!inherits(vs, "p2e_vs")) stop("vs must be a p2e_vs value set")
  mg <- ordinal_marginals(fit, newdata)
  a <- .vs_additive(vs)
  dec <- Reduce(`+`, lapply(seq_along(DIMS), function(d) drop(mg[[d]] %*% a$decrement[d, ])))
  u <- a$intercept - dec
  if (any(a$remainder != 0)) {
    rows <- seq_len(nrow(newdata))
    withr::with_seed(MC_SEED, .rng_kind = "Mersenne-Twister", .rng_normal_kind = "Inversion",
                     .rng_sample_kind = "Rejection", {
      for (b in split(rows, ceiling(rows / MVOP_MEAN_BLOCK))) {
        prof <- mvop_simulate_profiles(fit, newdata[b, , drop = FALSE], N_MC_MEAN)
        u[b] <- u[b] + rowMeans(matrix(a$remainder[.state_index(prof)], length(b)))
      }
    })
  }
  .clamp(u, vs)
}
