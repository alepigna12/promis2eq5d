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
