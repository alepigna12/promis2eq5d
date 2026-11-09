# Draw engine: an n x D utility matrix from a bundle.

# Rows per block (blocks hold whole keys); tests may lower it.
ENGINE_BLOCK <- 500L

# Row indices of each block. Keys are ordered by first appearance and added to
# the current block while it stays within `size` rows; a key larger than
# `size` forms its own block. Rows within a block are in ascending order.
#' @noRd
.block_rows <- function(key, size = ENGINE_BLOCK) {
  ord <- match(key, unique(key))
  counts <- tabulate(ord)
  block_of_key <- integer(length(counts))
  b <- 1L; used <- 0L
  for (k in seq_along(counts)) {
    if (used > 0L && used + counts[k] > size) { b <- b + 1L; used <- 0L }
    block_of_key[k] <- b
    used <- used + counts[k]
  }
  unname(split(seq_along(key), block_of_key[ord]))
}

.as_published_warning <- function(method) {
  why <- switch(method,
    D12_klapproth = paste("clamping the unbounded draws moves their mean below the point prediction,",
                          "and cohort intervals under-cover"),
    D03_clad = "about a quarter of the draws fall in the impossible gap below full health")
  warning("draw_type = \"as_published\" returns draws as published, which are biased: ", why,
          ". See ?predict_eq5d.", call. = FALSE)
}

.draw_matrix <- function(bundle, newdata, vs, D, key, draw_type = "corrected",
                         crosswalk = NULL, fits = bundle$boot) {
  if (!is.character(draw_type) || length(draw_type) != 1L || !draw_type %in% c("corrected", "as_published")) {
    stop("draw_type must be \"corrected\" or \"as_published\"", call. = FALSE)
  }
  m <- .p2e_method(bundle$method)
  n <- nrow(newdata)
  if (length(key) != n || !is.numeric(key) || anyNA(key) || any(key < 1) || any(key != round(key))) {
    stop("key must be a numeric vector with one value per row of newdata, no missing values, ",
         "and integer values of at least 1", call. = FALSE)
  }
  key <- as.integer(key)
  if (length(fits) < 1L) stop("fits must hold at least one fit", call. = FALSE)
  D <- as.integer(D)
  direct <- identical(m$family, "direct")
  if (direct) {
    if (!inherits(bundle$vs, "p2e_vs")) stop("direct bundle has no value set", call. = FALSE)
    if (!inherits(vs, "p2e_vs")) stop("a direct method needs the target vs to be a single value set", call. = FALSE)
    if (!identical(vs$code, bundle$vs$code) && is.null(crosswalk)) {
      stop("the bundle was fitted on value set '", bundle$vs$code, "' but the target is '", vs$code,
           "'; a crosswalk is needed", call. = FALSE)
    }
    # Still warns with a crosswalk: section 6 keeps g(mu-hat) for as_published, so the clamp and gap biases apply.
    if (draw_type == "as_published" && m$id %in% c("D12_klapproth", "D03_clad")) {
      .as_published_warning(m$id)
    }
    lo <- .vs_min(bundle$vs)
  }
  B <- length(fits)
  fit_index <- rep_len(seq_len(B), D)
  cols_of <- lapply(seq_len(B), function(b) which(fit_index == b))

  u <- matrix(NA_real_, n, D)
  profiles <- if (direct) NULL else array(NA_integer_, c(n, D, 5L))

  for (rows in .block_rows(key)) {
    nd_b <- newdata[rows, , drop = FALSE]
    key_b <- match(key[rows], unique(key[rows]))
    K_b <- max(key_b)
    nb <- length(rows)
    if (direct) {
      ub <- matrix(NA_real_, nb, D)
      for (b in seq_len(B)) {
        cols <- cols_of[[b]]
        if (!length(cols)) next
        ub[, cols] <- matrix(m$draw(fits[[b]], nd_b, length(cols), key_b, draw_type, lo), nb, length(cols))
      }
      ub <- .clamp(ub, bundle$vs)
      if (!is.null(crosswalk)) ub <- .apply_crosswalk(ub, crosswalk, vs)
      if (m$id == "D03_clad" && draw_type == "corrected") {
        V <- matrix(stats::runif(K_b * D), K_b, D)[key_b, , drop = FALSE]
        ub <- .round_gap(ub, vs, V)
      }
      u[rows, ] <- .clamp(ub, vs)
    } else {
      pb <- array(NA_integer_, c(nb, D, 5L))
      for (b in seq_len(B)) {
        cols <- cols_of[[b]]
        if (!length(cols)) next
        pb[, cols, ] <- m$simulate(fits[[b]], nd_b, length(cols), key_b)
      }
      profiles[rows, , ] <- pb
      u[rows, ] <- matrix(.score(pb, if (inherits(vs, "p2e_vs") || length(vs) == 1L) vs else vs[rows]), nb, D)
    }
  }
  list(u = u, profiles = profiles, fit_index = fit_index)
}
