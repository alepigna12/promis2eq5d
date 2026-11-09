# Draw corrections applied by the prediction engine after any crosswalk.

# Gap rounding: a draw u in the impossible gap b < u < 1, with b = 1 - .vs_gap(vs)
# (the highest attainable utility below full health is b), is moved to 1 with
# probability (u - b)/(1 - b) and to b otherwise, which keeps the mean. V is a
# matrix of uniforms the same size as u; other entries are unchanged.
.round_gap <- function(u, vs, V) {
  b <- 1 - .vs_gap(vs)
  inside <- u > b & u < 1
  out <- u
  out[inside] <- ifelse(V[inside] < (u[inside] - b) / (1 - b), 1, b)
  out
}

# Expected value of min(max(X, lo), hi) for X ~ N(m, sigma^2), in closed form,
# with a = (lo - m)/sigma and c = (hi - m)/sigma:
#   lo Phi(a) + hi (1 - Phi(c)) + m (Phi(c) - Phi(a)) + sigma (phi(a) - phi(c)).
# Vectorised over m.
.censored_normal_mean <- function(m, sigma, lo, hi) {
  a <- (lo - m) / sigma
  cc <- (hi - m) / sigma
  lo * stats::pnorm(a) + hi * (1 - stats::pnorm(cc)) +
    m * (stats::pnorm(cc) - stats::pnorm(a)) + sigma * (stats::dnorm(a) - stats::dnorm(cc))
}

# Location m of a normal N(m, sigma^2) censored to [lo, hi] whose censored
# mean E[min(max(m + sigma Z, lo), hi)] equals target. Vectorised over target;
# sigma is a scalar or as long as target. A target at or below lo gives -Inf
# and one at or above hi gives Inf: the draws are then point masses at the
# bound. Otherwise the root is found by uniroot on [lo - 10 sigma, hi + 10 sigma].
.censored_normal_location <- function(target, sigma, lo, hi) {
  if (anyNA(target) || anyNA(sigma)) stop("target and sigma must not be missing")
  n <- length(target)
  sigma <- rep_len(sigma, n)
  out <- numeric(n)
  out[target <= lo] <- -Inf
  out[target >= hi] <- Inf
  for (i in which(target > lo & target < hi)) {
    tg <- target[i]; s <- sigma[i]
    f <- function(m) .censored_normal_mean(m, s, lo, hi) - tg
    out[i] <- stats::uniroot(f, c(lo - 10 * s, hi + 10 * s), tol = 1e-10)$root
  }
  out
}
