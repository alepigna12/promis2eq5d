# Ported from iHEA_2027 R/functions/value_sets.R at freeze-ihea2027-abstract-v3 (a79e102); changes: country codes replaced by p2e_vs objects; the knots are fitted on the two utility vectors (no states table, no session cache, no from/to column names); apply_crosswalk takes explicit knots and has no identity or per-row-target paths; split the fit into .fit_knots(); added .choose_crosswalk_source()

# Weighted pool-adjacent-violators: non-decreasing least-squares fit to y
# (already ordered by x) with weights w > 0
#' @noRd
.pava <- function(y, w) {
  n <- length(y)
  val <- numeric(n); wt <- numeric(n); len <- integer(n); k <- 0L
  for (i in seq_len(n)) {
    k <- k + 1L; val[k] <- y[i]; wt[k] <- w[i]; len[k] <- 1L
    while (k > 1L && val[k - 1L] > val[k]) {
      nw <- wt[k - 1L] + wt[k]
      val[k - 1L] <- (wt[k - 1L] * val[k - 1L] + wt[k] * val[k]) / nw
      wt[k - 1L] <- nw; len[k - 1L] <- len[k - 1L] + len[k]
      k <- k - 1L
    }
  }
  rep(val[seq_len(k)], len[seq_len(k)])
}

# Knots data.frame(x, y) from paired vectors x, y with weights w. Rows with
# weight 0 are dropped, tied x are pooled at their weighted mean y, then a
# weighted PAVA is run. xr and yr are c(min, max) of the full from / to ranges:
# the fitted y is clamped into yr (rounding can push a pooled mean an ulp
# outside it) and the knots are anchored at (xr[1], yr[1]) and (xr[2], yr[2])
# if missing.
#' @noRd
.fit_knots <- function(x, y, w, xr, yr) {
  keep <- w > 0
  x <- x[keep]; y <- y[keep]; w <- w[keep]
  ux <- sort(unique(x))
  s  <- rowsum(cbind(w, w * y), match(x, ux))
  k <- data.frame(x = ux, y = .pava(s[, 2] / s[, 1], s[, 1]))
  k$y <- pmin(pmax(k$y, yr[1]), yr[2])
  if (k$x[1] > xr[1]) k <- rbind(data.frame(x = xr[1], y = yr[1]), k)
  if (k$x[nrow(k)] < xr[2]) k <- rbind(k, data.frame(x = xr[2], y = yr[2]))
  rownames(k) <- NULL
  k
}

# weights: NULL (equal) or a finite, non-negative numeric vector of 3,125
# state weights (state order), not all zero
#' @noRd
.check_weights <- function(weights) {
  if (is.null(weights)) return(rep(1, 3125))
  if (!is.numeric(weights) || length(weights) != 3125L || !all(is.finite(weights)) ||
      any(weights < 0)) {
    stop("weights must be finite and non-negative numbers, one per state (3125)")
  }
  if (!any(weights > 0)) stop("weights must not all be zero")
  weights
}

# Isotonic crosswalk knots data.frame(x, y): the fitted `to` utility at each
# distinct `from` utility, over the 3,125 states. Attributes "from" and "to"
# carry the value-set codes.
#' @noRd
.fit_crosswalk <- function(from_vs, to_vs, weights = NULL) {
  from_vs <- p2e_value_set(from_vs); to_vs <- p2e_value_set(to_vs)
  w <- .check_weights(weights)
  k <- .fit_knots(from_vs$u, to_vs$u, w, range(from_vs$u), range(to_vs$u))
  attr(k, "from") <- from_vs$code; attr(k, "to") <- to_vs$code
  k
}

# Map utilities u (vector or matrix) on the from value set to to_vs: linear
# interpolation over the knots, constant beyond the end knots, then clamped to
# [.vs_min(to_vs), 1].
#' @noRd
.apply_crosswalk <- function(u, knots, to_vs) {
  to_vs <- p2e_value_set(to_vs)
  if (!is.data.frame(knots) || !identical(names(knots)[1:2], c("x", "y")) ||
      anyNA(knots$x) || anyNA(knots$y) || nrow(knots) < 2 ||
      any(diff(knots$x) <= 0) || any(diff(knots$y) < 0)) {
    stop("knots must be data.frame(x, y) with x strictly increasing and y non-decreasing")
  }
  if (anyNA(u)) stop("u contains missing values")
  out <- u
  out[] <- stats::approx(knots$x, knots$y, xout = as.vector(u), rule = 2, ties = "ordered")$y
  .clamp(out, to_vs)
}

# The source value set (of `sources`) whose crosswalk to to_vs has the lowest
# weighted RMSE over the 3,125 states; the first in `sources` order on a tie.
#' @noRd
.choose_crosswalk_source <- function(to_vs, sources = c("FR", "DE", "UK"), weights = NULL) {
  to_vs <- p2e_value_set(to_vs)
  w <- .check_weights(weights)
  nw <- w / sum(w)
  rmse <- vapply(sources, function(s) {
    f <- p2e_value_set(s)
    mapped <- .apply_crosswalk(f$u, .fit_crosswalk(f, to_vs, w), to_vs)
    sqrt(sum(nw * (mapped - to_vs$u)^2))
  }, numeric(1))
  sources[which.min(rmse)]
}
