# The p2e_result class: print, summary and as.matrix methods.

# The result as a plain data frame
.plain_result <- function(x) {
  attr(x, "p2e_info") <- NULL
  class(x) <- "data.frame"
  x
}

.result_info <- function(x) {
  info <- attr(x, "p2e_info")
  if (is.null(info)) stop("not a p2e_result: the p2e_info attribute is missing.", call. = FALSE)
  info
}

.result_route <- function(info) {
  if (identical(info$family, "response")) return("response method, scored under the value set")
  if (!is.na(info$crosswalk)) return(paste0("iso crosswalk from ", info$crosswalk))
  if (isTRUE(info$retrain)) return("refitted on the value set")
  "shipped fit"
}

#' Print, summarise and convert prediction results
#'
#' Methods for the result of [predict_eq5d()].
#'
#' @param x,object A `p2e_result`.
#' @param n Number of rows to print.
#' @param level `summary`, individual mode: the level of the draw interval.
#'   Cohort summaries use the 95% level.
#' @param ... Not used.
#'
#' @return `print` returns `x` invisibly. `summary` returns a data frame:
#'   mean mode, `n`, `mean`, `sd`, `min`, `max` of `utility`; individual mode,
#'   one row per input row with `mean`, `sd`, `lower` and `upper`; cohort mode,
#'   one row per group with `n`, `mean`, `lower` and `upper`. `as.matrix`
#'   returns the draws (rows or groups by draws) and is an error in mean mode.
#' @seealso [predict_eq5d()]
#' @name p2e_result
NULL

#' @rdname p2e_result
#' @export
print.p2e_result <- function(x, n = 6, ...) {
  info <- .result_info(x)
  vs <- info$value_set
  cat("<p2e_result> ", nrow(x), " rows\n", sep = "")
  cat("method: ", info$method, "\n", sep = "")
  cat("set: ", info$set, " (", SET_DESCRIPTIONS[[info$set]], ")\n", sep = "")
  cat("value set: ", vs$name, " (", vs$code, ")\n", sep = "")
  cat("route: ", .result_route(info), "\n", sep = "")
  cat("mode: ", info$mode, "\n", sep = "")
  cat("draws: ", if (is.na(info$draws)) "none (mean mode)" else info$draws, "\n", sep = "")
  df <- .plain_result(x)
  print(utils::head(df, n))
  if (nrow(df) > n) cat("... with ", nrow(df) - n, " more rows\n", sep = "")
  invisible(x)
}

# Group index (by first appearance of the group columns) of a cohort result
.cohort_index <- function(df) {
  gcols <- setdiff(names(df), c(".draw", "n", "mean_utility"))
  lab <- if (length(gcols)) do.call(paste, c(lapply(df[gcols], as.character), sep = "\r")) else
    rep("", nrow(df))
  key <- match(lab, unique(lab))
  first <- match(seq_len(max(key)), key)
  labels <- if (length(gcols)) do.call(paste, c(lapply(df[gcols], as.character), sep = ", ")) else
    "all"
  list(gcols = gcols, key = key, first = first, labels = labels[first])
}

#' @rdname p2e_result
#' @export
summary.p2e_result <- function(object, level = 0.90, ...) {
  info <- .result_info(object)
  df <- .plain_result(object)
  if (info$mode == "mean") {
    u <- df$utility
    return(data.frame(n = length(u), mean = mean(u), sd = stats::sd(u), min = min(u), max = max(u)))
  }
  if (!(is.numeric(level) && length(level) == 1L && level > 0 && level < 1)) {
    stop("level must be a single number strictly between 0 and 1.", call. = FALSE)
  }
  if (info$mode == "individual") {
    u <- as.matrix(object)
    a <- (1 - level) / 2
    first <- match(sort(unique(df$.row)), df$.row)
    idc <- setdiff(names(df), c(".row", ".draw", "utility", DIM_NAMES))
    out <- df[first, c(".row", idc), drop = FALSE]
    out$mean <- unname(rowMeans(u))
    out$sd <- unname(apply(u, 1, stats::sd))
    out$lower <- unname(apply(u, 1, stats::quantile, a))
    out$upper <- unname(apply(u, 1, stats::quantile, 1 - a))
    rownames(out) <- NULL
    return(out)
  }
  ix <- .cohort_index(df)
  a <- (1 - COHORT_LEVEL) / 2
  out <- df[ix$first, c(ix$gcols, "n"), drop = FALSE]
  by_g <- split(df$mean_utility, factor(ix$key, levels = seq_along(ix$first)))
  out$mean <- unname(vapply(by_g, mean, numeric(1)))
  out$lower <- unname(vapply(by_g, stats::quantile, numeric(1), a))
  out$upper <- unname(vapply(by_g, stats::quantile, numeric(1), 1 - a))
  rownames(out) <- NULL
  out
}

#' @rdname p2e_result
#' @export
as.matrix.p2e_result <- function(x, ...) {
  # a subset such as x[1:5, ] loses p2e_info: treat it as a plain data frame
  if (is.null(attr(x, "p2e_info"))) return(as.matrix(.plain_result(x), ...))
  info <- .result_info(x)
  if (info$mode == "mean") {
    stop("a result in mean mode has no draws; use mode = \"individual\" or \"cohort\".", call. = FALSE)
  }
  df <- .plain_result(x)
  D <- max(df$.draw)
  if (info$mode == "individual") {
    rows <- sort(unique(df$.row))
    m <- matrix(NA_real_, length(rows), D)
    m[cbind(match(df$.row, rows), df$.draw)] <- df$utility
    dimnames(m) <- list(NULL, NULL)
    return(m)
  }
  ix <- .cohort_index(df)
  m <- matrix(NA_real_, length(ix$first), D)
  m[cbind(ix$key, df$.draw)] <- df$mean_utility
  dimnames(m) <- list(ix$labels, NULL)
  m
}
