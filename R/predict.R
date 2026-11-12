# The main function: predict_eq5d().

DIM_NAMES <- c("MO", "SC", "UA", "PD", "AD")
SHIPPED_TARGETS <- c("FR", "DE", "UK")
SET_DESCRIPTIONS <- c(
  S1 = "items + pain intensity + age/sex",
  S2 = "items + pain intensity",
  S3 = "T-scores + pain intensity + age/sex",
  S4 = "T-scores + pain intensity",
  S3K = "T-scores + age/sex"
)

#' Predict EQ-5D-5L utilities from PROMIS-29 profiles
#'
#' Maps PROMIS-29 profiles to EQ-5D-5L utilities, with uncertainty, under the
#' value set of any country in the 'eq5d' package (or your own value set).
#' See [p2e_columns()] for the columns the function reads. Nothing is
#' imputed: missing values in the columns of the chosen predictor set are an
#' error.
#'
#' @section Family and method:
#' `family` is `"response"` (the profile is predicted and scored under the
#' value set; any country works) or `"direct"` (the utility is predicted on
#' the scale of one value set). `method` is one of `"R06_mvop_dim"` and
#' `"R09_chain"` (response) or `"D03_clad"` and `"D12_klapproth"` (direct);
#' short aliases (`"R06"`, `"mvop"`, `"chain"`, `"clad"`, `"klapproth"`) are
#' accepted. With neither argument the default response method
#' (`"R06_mvop_dim"`) is used; with `family` only, that family's default
#' (`"D03_clad"` for direct); with `method` only, the family follows from the
#' method; with both, they must agree.
#'
#' @section Predictor set:
#' Each method uses its preferred predictor set among those the data support
#' (items, T-scores, pain intensity and age/sex; see [p2e_columns()]). `set`
#' forces one of `"S1"` to `"S4"`. The chosen set is reported in a message and
#' in the result. The package ships no PROMIS T-score conversion table: T-score
#' sets need T-score columns or a `tscore_table` (see [p2e_tscores()]).
#'
#' @section Modes:
#' * `"mean"`: one expected utility per row (`utility`), with `interval = TRUE`
#'   a CV+ conformal interval (`lower`, `upper`) at `level`, and with
#'   `dimensions = TRUE` the most likely level of each dimension (`MO`, `SC`,
#'   `UA`, `PD`, `AD`; response methods only).
#' * `"individual"`: `draws` draws per row (`.row`, `.draw`, `utility`), with
#'   `dimensions = TRUE` the sampled profile.
#' * `"cohort"`: for each draw, the mean utility of each `group` (all rows when
#'   `group` is `NULL`), with the group size `n`.
#'
#' Draw d uses the same parameter draw (the same bootstrap fit) for every row,
#' group and time point in a call. This keeps draws aligned for QALY
#' calculations: use one call with `group`, or separate calls with the same
#' `draws`, which are aligned whatever the seed. Give separate calls different
#' seeds (or none): with the same `seed`, rows in the same position of two
#' calls share their individual random numbers, which understates the
#' uncertainty of differences between them. [summary.p2e_result()] and [as.matrix.p2e_result()] work on
#' the draws.
#'
#' @section Reaching your country:
#' Response methods score their predicted profiles under any value set;
#' `retrain` has no effect for them. Direct methods ship fits for `"FR"`,
#' `"DE"` and `"UK"`. For any other country, with `retrain = FALSE`, the
#' prediction and each draw are mapped to the country by an isotonic crosswalk
#' from one shipped value set (by default the one with the lowest weighted
#' RMSE; `crosswalk_from` overrides). With `retrain = TRUE` the method is
#' refitted on the shipped training data rescored under the country's value
#' set (slow; `boot` sets the number of bootstrap refits). Interval coverage
#' for a value set never seen in training relies on exchangeability.
#'
#' @section D03 mean or median:
#' With `statistic = "mean"` (default) the D03 utility is the mean of its
#' corrected residual draws (a fixed Monte Carlo seed makes it reproducible),
#' so mean mode agrees with cohort means. `statistic = "median"` gives the
#' censored median, which is the pipeline's D03 point prediction. Conformal
#' intervals for D03 stay centred on the median. `statistic` has no effect on
#' other methods.
#'
#' @section Draws for direct methods:
#' `draw_type = "corrected"` (default) draws with the mean preserved:
#' D12 from a censored normal, D03 with mean-preserving rounding of draws in
#' the gap below full health (this does not fully restore the full-health
#' share; use a response method when it matters). `"as_published"` returns the
#' draws as published, which are biased, and warns once per call. In mean mode
#' it returns the pipeline mean, the clamped (crosswalked) point prediction.
#'
#' @section Several rows per person:
#' With `id` naming a person column and `within_id = "independent"` (default)
#' the person-level random part of each draw is independent on each row. With
#' `within_id = "common"` rows with the same `id` share their random numbers
#' within each draw, so a person's errors are perfectly rank-correlated across
#' rows. The training data are cross-sectional, so the true correlation cannot
#' be estimated: treat `"common"` as a sensitivity analysis. Parameter draws
#' stay aligned across all rows in both settings.
#'
#' @section Trademarks:
#' EQ-5D\Sexpr[results=rd]{"\\\\enc{\u2122}{(TM)}"} is a trade mark of the EuroQol Research Foundation.
#'
#' PROMIS\Sexpr[results=rd]{"\\\\enc{\u00ae}{(R)}"} is a registered trademark of the U.S. Department of Health and Human Services. This package is not affiliated with or endorsed by HHS, the PROMIS Health Organization, HealthMeasures or the EuroQol Research Foundation.
#'
#' @param data A data frame (or tibble) with one row per person or per person
#'   and time point.
#' @param country A country name or code known to 'eq5d' (see
#'   [p2e_value_set()]), a `p2e_vs` object, or a table of 3,125 states and
#'   utilities (passed to [p2e_value_set()]).
#' @param family `NULL`, `"response"` or `"direct"`.
#' @param method `NULL`, or a method id or alias.
#' @param mode `"mean"`, `"individual"` or `"cohort"`.
#' @param group Cohort mode: names of the columns that define the groups.
#'   Groups are ordered by first appearance. Missing values are an error.
#' @param draws Number of draws in the draw modes.
#' @param interval Mean mode: add the CV+ conformal interval.
#' @param level Level of the conformal interval.
#' @param dimensions Also return the dimension levels (response methods only).
#' @param set `NULL` or a predictor set, `"S1"` to `"S4"`.
#' @param retrain Direct methods and a country without a shipped fit: refit
#'   instead of using the crosswalk.
#' @param draw_type `"corrected"` or `"as_published"` (D03 and D12 only).
#' @param id Name of a column identifying a person across rows. Missing values
#'   are an error with `within_id = "common"`.
#' @param within_id `"independent"` or `"common"`.
#' @param statistic Mean mode, D03 only: `"mean"` or `"median"`.
#' @param crosswalk_from Direct methods with a crosswalk: the source value set,
#'   `"FR"`, `"DE"` or `"UK"`.
#' @param boot `retrain = TRUE`: number of bootstrap refits (default all).
#' @param tscore_table A user-supplied raw-to-T-score table (see
#'   [p2e_tscores()]).
#' @param seed `NULL` or a number that makes the draws reproducible.
#'
#' @return A data frame of class `p2e_result` with the attribute `p2e_info`
#'   (method, family, set, value set, mode, draws, `B` (the number of bootstrap
#'   fits used), `boot_attempted` and `boot_failed` (the bootstrap fits
#'   attempted when the models were fitted, and those that failed or did not
#'   converge; intervals and draws are conditional on convergence),
#'   crosswalk source, retrain, draw type, `within_id`, seed and package
#'   version). In the draw modes `print()` reports the bootstrap counts and
#'   `summary()` carries them as the attributes `boot_attempted` and
#'   `boot_failed`.
#' @seealso [p2e_columns()], [p2e_value_set()], [p2e_tscores()]
#' @examples
#' \donttest{
#' pro <- data.frame(age = 54, female = 1, Global07 = 4)
#' pro[unname(subset(p2e_columns(), type == "item")$column)] <- 2L
#' predict_eq5d(pro, "UK")
#' predict_eq5d(pro, "Japan", mode = "individual", draws = 200, seed = 1)
#' }
#' @export
predict_eq5d <- function(data, country, family = NULL, method = NULL,
                         mode = c("mean", "individual", "cohort"), group = NULL,
                         draws = 1000, interval = FALSE, level = 0.90,
                         dimensions = FALSE, set = NULL, retrain = FALSE,
                         draw_type = c("corrected", "as_published"), id = NULL,
                         within_id = c("independent", "common"),
                         statistic = c("mean", "median"), crosswalk_from = NULL,
                         boot = NULL, tscore_table = NULL, seed = NULL) {
  mode <- match.arg(mode)
  draw_type <- match.arg(draw_type)
  within_id <- match.arg(within_id)
  statistic <- match.arg(statistic)
  .check_flag(interval, "interval"); .check_flag(dimensions, "dimensions")
  .check_flag(retrain, "retrain")
  if (!is.null(seed) && !(is.numeric(seed) && length(seed) == 1L && is.finite(seed))) {
    stop("seed must be NULL or a single number.", call. = FALSE)
  }
  draw_mode <- mode != "mean"
  if (draw_mode && !(is.numeric(draws) && length(draws) == 1L && is.finite(draws) && draws >= 1)) {
    stop("draws must be a single number of at least 1.", call. = FALSE)
  }
  D <- as.integer(draws)
  if (interval && draw_mode) {
    warning("interval is only used in mean mode; it is ignored.", call. = FALSE)
    interval <- FALSE
  }
  if (interval && !(is.numeric(level) && length(level) == 1L && !is.na(level) &&
                    level > 0 && level < 1)) {
    stop("level must be a single number strictly between 0 and 1.", call. = FALSE)
  }
  if (!is.null(group) && mode != "cohort") {
    warning("group is only used in cohort mode; it is ignored.", call. = FALSE)
  }

  data <- .plain_frame(data)
  if (nrow(data) == 0L) stop("`data` has no rows.", call. = FALSE)
  n <- nrow(data)

  vs <- p2e_value_set(country)
  m <- .resolve_method(family, method)
  if (dimensions && identical(m$family, "direct")) {
    stop("dimension levels need a response method (R06_mvop_dim or R09_chain).", call. = FALSE)
  }
  .require_pkgs(m$predict_needs, paste0("predict with ", m$id))
  direct <- identical(m$family, "direct")

  idcol <- .check_key_col(data, id, "id", missing_error = within_id == "common")
  gcols <- if (mode == "cohort") .check_group_cols(data, group) else NULL
  key <- if (!is.null(id) && within_id == "common") match(data[[id]], unique(data[[id]])) else seq_len(n)

  std <- .standardise_input(data, tscore_table)
  set <- tryCatch(.choose_set(m$id, std, set), error = function(e) {
    # incomplete items and no T-scores: say why the table did not help
    miss <- setdiff(ITEMS, names(std))
    if (!is.null(tscore_table) && length(miss) && !all(TSCORES %in% names(std))) {
      stop("a T-score table was given, but the data have neither all 28 items nor all T-score ",
           "columns; missing item column(s): ",
           paste(vapply(miss, .user_name, character(1)), collapse = ", "), ".", call. = FALSE)
    }
    stop(e)
  })
  message("Using predictor set ", set, " (", SET_DESCRIPTIONS[[set]], ")")

  if (!direct && retrain) message("retrain has no effect for response methods")
  route <- .resolve_bundle(m, set, vs, direct, retrain, crosswalk_from, boot, tscore_table,
                           interval, mode)
  bundle <- route$bundle; knots <- route$knots

  run <- function(expr) if (is.null(seed)) expr else withr::with_seed(seed, expr)

  if (mode == "mean") {
    u <- .mean_utility(m, bundle, std, vs, knots, statistic, draw_type)
    out <- data.frame(.row = seq_len(n))
    if (!is.null(idcol)) out[[id]] <- data[[id]]
    out$utility <- u
    if (interval) {
      ci <- .conformal(bundle, std, vs, level, crosswalk = knots)
      out$lower <- ci$lower; out$upper <- ci$upper
    }
    if (dimensions) {
      mg <- m$marginals(bundle$point, std)
      for (j in seq_along(DIMS)) out[[DIM_NAMES[j]]] <- max.col(mg[[DIMS[j]]], ties.method = "first")
    }
  } else {
    res <- run(.draw_matrix(bundle, std, vs, D, key, draw_type, crosswalk = knots))
    if (mode == "individual") {
      out <- data.frame(.row = rep(seq_len(n), D))
      if (!is.null(idcol)) out[[id]] <- data[[id]][out$.row]
      out$.draw <- rep(seq_len(D), each = n)
      out$utility <- as.vector(res$u)
      if (dimensions) {
        for (j in seq_along(DIMS)) out[[DIM_NAMES[j]]] <- as.vector(res$profiles[, , j])
      }
    } else {
      out <- .cohort_draws(data, gcols, res$u)
    }
  }

  info <- list(
    method = m$id, family = m$family, set = set,
    value_set = list(code = vs$code, name = vs$name, version = vs$version, source = vs$source),
    mode = mode, draws = if (draw_mode) D else NA_integer_, B = length(bundle$boot),
    boot_attempted = if (is.null(bundle$boot_attempted)) NA_integer_ else as.integer(bundle$boot_attempted),
    boot_failed = if (is.null(bundle$boot_failed)) NA_integer_ else as.integer(bundle$boot_failed),
    crosswalk = if (is.null(knots)) NA_character_ else route$source,
    retrain = route$retrained, draw_type = draw_type, within_id = within_id,
    seed = if (is.null(seed)) NA else seed,
    version = utils::packageVersion("promis2eq5d")
  )
  structure(out, p2e_info = info,
            p2e_cols = list(id = idcol, group = gcols),
            class = c("p2e_result", "data.frame"))
}

.check_flag <- function(x, name) {
  if (!(is.logical(x) && length(x) == 1L && !is.na(x))) {
    stop(name, " must be TRUE or FALSE.", call. = FALSE)
  }
}

# A plain data frame whatever the input class (data.frame, tibble, grouped
# tibble); no tibble or dplyr dependency.
.plain_frame <- function(data) {
  if (!is.data.frame(data)) stop("`data` must be a data frame.", call. = FALSE)
  cols <- unclass(data)
  nr <- nrow(data)
  attributes(cols) <- list(names = names(data))
  out <- as.data.frame(cols, stringsAsFactors = FALSE, optional = TRUE)
  if (!length(cols)) out <- data.frame(row.names = seq_len(nr))
  out
}

# Family/method resolution (the 0/1/2 table)
.resolve_method <- function(family, method) {
  if (!is.null(family)) {
    if (!(is.character(family) && length(family) == 1L && family %in% c("response", "direct"))) {
      stop("family must be \"response\" or \"direct\".", call. = FALSE)
    }
  }
  if (is.null(method)) {
    return(.p2e_method(switch(if (is.null(family)) "response" else family,
                              response = "R06_mvop_dim", direct = "D03_clad")))
  }
  m <- .p2e_method(method)
  if (!is.null(family) && !identical(family, m$family)) {
    stop(m$id, " is a ", m$family, " method, not a ", family, " method.", call. = FALSE)
  }
  m
}

# Validate an id or group column: it must exist; missing values are an error
# when missing_error is TRUE. Returns the column name or NULL.
.check_key_col <- function(data, col, arg, missing_error = TRUE) {
  if (is.null(col)) return(NULL)
  if (!(is.character(col) && length(col) == 1L && !is.na(col))) {
    stop(arg, " must be the name of one column.", call. = FALSE)
  }
  if (!col %in% names(data)) stop(arg, " column '", col, "' is not in `data`.", call. = FALSE)
  na <- which(is.na(data[[col]]))
  if (length(na) && missing_error) {
    stop(arg, " column '", col, "' has missing values at ", .row_text(na), ".", call. = FALSE)
  }
  col
}

.check_group_cols <- function(data, group) {
  if (is.null(group)) return(NULL)
  if (!is.character(group) || !length(group) || anyNA(group)) {
    stop("group must be a character vector of column names.", call. = FALSE)
  }
  for (g in group) .check_key_col(data, g, "group")
  group
}

# The fit bundle, the crosswalk knots (or NULL) and the route used
.resolve_bundle <- function(m, set, vs, direct, retrain, crosswalk_from, boot, tscore_table,
                            interval, mode) {
  if (!direct) {
    return(list(bundle = .load_bundle(m$id, set, NA), knots = NULL, source = NA_character_,
                retrained = FALSE))
  }
  shipped <- identical(vs$source, "eq5d") && vs$code %in% SHIPPED_TARGETS
  if (shipped) {
    if (retrain) message("retrain has no effect: ", vs$code, " has a shipped fit")
    return(list(bundle = .load_bundle(m$id, set, vs$code), knots = NULL, source = NA_character_,
                retrained = FALSE))
  }
  if (retrain) {
    b <- .retrain_bundle(m$id, set, vs, boot = boot, tscore_table = tscore_table,
                         interval = interval, mode = mode)
    return(list(bundle = b, knots = NULL, source = NA_character_, retrained = TRUE))
  }
  uk <- .load_bundle(m$id, set, "UK")
  freq <- uk$state_freq
  src <- if (is.null(crosswalk_from)) {
    .choose_crosswalk_source(vs, weights = freq)
  } else {
    if (!(is.character(crosswalk_from) && length(crosswalk_from) == 1L &&
          crosswalk_from %in% SHIPPED_TARGETS)) {
      stop("crosswalk_from must be one of the shipped value sets: ",
           paste(SHIPPED_TARGETS, collapse = ", "), ".", call. = FALSE)
    }
    crosswalk_from
  }
  b <- if (identical(src, "UK")) uk else .load_bundle(m$id, set, src)
  knots <- .fit_crosswalk(b$vs, vs, freq)
  message("Using the iso crosswalk from ", src, " to ", vs$name,
          " (direct method; set retrain = TRUE to refit)")
  list(bundle = b, knots = knots, source = src, retrained = FALSE)
}

# One expected utility per row of newdata (std)
.mean_utility <- function(m, bundle, newdata, vs, knots, statistic, draw_type) {
  if (identical(m$family, "response")) return(m$mean_utility(bundle$point, newdata, vs))
  n <- nrow(newdata)
  g <- function(p) {
    p <- .clamp(p, bundle$vs)
    if (!is.null(knots)) p <- .apply_crosswalk(p, knots, vs)
    .clamp(p, vs)
  }
  p <- m$point(bundle$point, newdata)
  if (draw_type == "as_published") {
    .as_published_warning(m$id)
    return(g(p))
  }
  if (m$id == "D12_klapproth" && is.null(knots)) return(.clamp(p, bundle$vs))
  if (m$id == "D03_clad" && statistic == "median") return(g(p))
  # row blocks bound the memory (n x N_MC_MEAN per block); one stream
  rows <- seq_len(n)
  out <- numeric(n)
  .with_mc_seed({
    for (b in split(rows, ceiling(rows / MC_MEAN_BLOCK))) {
      out[b] <- rowMeans(.draw_matrix(bundle, newdata[b, , drop = FALSE], vs, N_MC_MEAN,
                                      seq_along(b), "corrected", crosswalk = knots,
                                      fits = list(bundle$point))$u)
    }
  })
  out
}

# Rows per block of the D03/D12 Monte Carlo mean (tests may lower it)
MC_MEAN_BLOCK <- 500L

# Cohort draws: group means of the draws matrix u (n x D)
.cohort_draws <- function(data, gcols, u) {
  n <- nrow(u); D <- ncol(u)
  gkey <- if (is.null(gcols)) rep(1L, n) else {
    lab <- do.call(paste, c(lapply(data[gcols], as.character), sep = "\r"))
    match(lab, unique(lab))
  }
  G <- max(gkey)
  first <- match(seq_len(G), gkey)
  means <- t(vapply(seq_len(G), function(g) colMeans(u[gkey == g, , drop = FALSE]), numeric(D)))
  dim(means) <- c(G, D)
  out <- if (is.null(gcols)) data.frame(row.names = seq_len(G * D)) else {
    cols <- lapply(data[gcols], function(x) drop_unused_levels(x)[first][rep(seq_len(G), each = D)])
    as.data.frame(cols, stringsAsFactors = FALSE, optional = TRUE)
  }
  out$.draw <- rep(seq_len(D), times = G)
  out$n <- rep(tabulate(gkey, G), each = D)
  out$mean_utility <- as.vector(t(means))
  out
}

# Factor columns lose the levels that no row uses
drop_unused_levels <- function(x) if (is.factor(x)) droplevels(x) else x
