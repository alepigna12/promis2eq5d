# Input layer: user-facing columns, validation and predictor-set choice.

#' Input columns accepted by the package
#'
#' Lists the data columns that the prediction functions read, the type and
#' allowed range of each, and the predictor sets that need each column. Item
#' columns use the HealthMeasures item IDs. Pain intensity may also be given
#' as `P9`, and age as `AGE`. Sex is given either as `female` (0/1) or as
#' `sex` (for example `"F"`/`"M"` or `"female"`/`"male"`; `"woman"`/`"man"` are also accepted, case-insensitive).
#' Missing values are allowed in the data; they are an error only in a column
#' of the predictor set that is used.
#'
#' @return A data frame with columns `column` (the name to use in your data),
#'   `internal` (the name used inside the package), `type` (`"item"`,
#'   `"pain"`, `"tscore"` or `"covariate"`), `range` and `required_by` (the
#'   predictor sets that use the column).
#' @export
#' @examples
#' p2e_columns()
#' subset(p2e_columns(), type == "covariate")
p2e_columns <- function() {
  req <- function(internal) {
    vapply(internal, function(v) {
      paste(names(PREDICTOR_SETS)[vapply(PREDICTOR_SETS, function(s) v %in% s,
                                         logical(1))], collapse = ", ")
    }, character(1), USE.NAMES = FALSE)
  }
  mk <- function(column, internal, type, range) {
    data.frame(column = column, internal = internal, type = type, range = range,
               required_by = req(internal), stringsAsFactors = FALSE)
  }
  rbind(
    mk(unname(ITEM_IDS), names(ITEM_IDS), "item", "1-5"),
    mk(PAIN_ID, "P9", "pain", "0-10"),
    mk(TSCORES, TSCORES, "tscore", "numeric"),
    mk(c("age", "female", "sex"), c("AGE", "female", "female"), "covariate",
       c("18-100 (warning outside)", "0 or 1", "F/M, female/male"))
  )
}

# Validate one column and return it as numeric.
.check_numeric <- function(x, col) {
  if (is.factor(x)) {
    stop("column '", col, "' is a factor; supply a numeric column.", call. = FALSE)
  }
  if (!is.numeric(x) && !is.logical(x)) {
    stop("column '", col, "' must be numeric (not ", class(x)[1], ").",
         call. = FALSE)
  }
  x <- as.numeric(x)
  inf <- which(is.infinite(x))
  if (length(inf)) {
    stop("column '", col, "' has non-finite values at ", .row_text(inf), ".",
         call. = FALSE)
  }
  x
}

.row_text <- function(i) {
  shown <- paste(utils::head(i, 5), collapse = ", ")
  if (length(i) > 5) shown <- paste0(shown, ", ...")
  paste0(if (length(i) == 1) "row " else "rows ", shown)
}

.check_range <- function(x, col, lo, hi, integer = FALSE) {
  bad <- which(!is.na(x) & (x < lo | x > hi | (integer & x != round(x))))
  if (length(bad)) {
    stop("column '", col, "' must be ", if (integer) "an integer " else "",
         "in ", lo, "-", hi, "; offending values at ", .row_text(bad), ".",
         call. = FALSE)
  }
  x
}

# Pick the column for one input with a user-facing name and an alias.
.pick <- function(data, name, alias) {
  have <- intersect(c(name, alias), names(data))
  if (length(have) > 1) {
    stop("give only one of the columns ",
         paste0("'", have, "'", collapse = " and "), ".", call. = FALSE)
  }
  if (length(have)) have else NA_character_
}

.sex_to_female <- function(x, col) {
  if (is.factor(x)) x <- as.character(x)
  if (is.character(x)) {
    k <- tolower(trimws(x))
    out <- ifelse(k %in% c("f", "female", "woman"), 1,
                  ifelse(k %in% c("m", "male", "man"), 0, NA_real_))
    bad <- which(!is.na(x) & is.na(out))
    if (length(bad)) {
      stop("column '", col, "' has values that are not F/M or female/male at ",
           .row_text(bad), "; or supply `female` as 0/1.", call. = FALSE)
    }
    return(out)
  }
  stop("column '", col, "' must be text such as \"F\"/\"M\" or ",
       "\"female\"/\"male\"; or supply `female` as 0/1.", call. = FALSE)
}

# Map user columns to the internal names and validate them. Returns a data
# frame with the recognised columns only, in the row order of `data`.
.standardise_input <- function(data, tscore_table = NULL) {
  if (!is.data.frame(data)) stop("`data` must be a data frame.", call. = FALSE)
  data <- as.data.frame(data, stringsAsFactors = FALSE)
  out <- list()

  for (nm in ITEMS) {
    id <- ITEM_IDS[[nm]]
    if (!id %in% names(data)) next
    out[[nm]] <- .check_range(.check_numeric(data[[id]], id), id, 1, 5,
                              integer = TRUE)
  }

  src <- .pick(data, PAIN_ID, "P9")
  if (!is.na(src)) {
    out$P9 <- .check_range(.check_numeric(data[[src]], src), src, 0, 10)
  }

  for (nm in TSCORES) {
    if (nm %in% names(data)) out[[nm]] <- .check_numeric(data[[nm]], nm)
  }

  src <- .pick(data, "age", "AGE")
  if (!is.na(src)) {
    a <- .check_numeric(data[[src]], src)
    td <- tryCatch(.training_data(), error = function(e) NULL)
    rng <- if (is.null(td)) c(18, 100) else c(floor(min(td$AGE)), ceiling(max(td$AGE)))
    odd <- which(!is.na(a) & (a < rng[1] | a > rng[2]))
    if (length(odd)) {
      warning("column '", src, "' has values outside ", rng[1], "-", rng[2], " at ", .row_text(odd),
              "; predictions outside the training age range are extrapolations.",
              call. = FALSE)
    }
    out$AGE <- a
  }

  if (all(c("female", "sex") %in% names(data))) {
    stop("give only one of the columns 'female' and 'sex'.", call. = FALSE)
  }
  fem <- intersect(c("female", "sex"), names(data))
  if (length(fem)) {
    out$female <- if (fem == "female") {
      .check_range(.check_numeric(data[[fem]], fem), fem, 0, 1, integer = TRUE)
    } else {
      .sex_to_female(data[[fem]], fem)
    }
  }

  std <- data.frame(row.names = seq_len(nrow(data)))
  for (nm in names(out)) std[[nm]] <- out[[nm]]
  rownames(std) <- rownames(data)

  # The table scores the user's data only when all 28 items are present. With
  # incomplete items, T_* columns already in the data are used as given (the
  # table is then used only to score the training data when retraining).
  miss <- setdiff(ITEMS, names(std))
  if (!is.null(tscore_table) && !length(miss)) {
    ts <- .tscores_internal(std, tscore_table)
    replaced <- intersect(TSCORES, names(std))
    if (length(replaced)) {
      warning("replacing the supplied column(s) ", paste(replaced, collapse = ", "),
              " with T-scores computed from `tscore_table`.", call. = FALSE)
    }
    for (nm in TSCORES) std[[nm]] <- ts[[nm]]
  }
  std
}

# Sets whose columns are all present in `std`.
.available_sets <- function(std) {
  names(PREDICTOR_SETS)[vapply(PREDICTOR_SETS, function(s) all(s %in% names(std)),
                               logical(1))]
}

# Choose the predictor set for a method: the first available set in the
# method's preferred order, or the requested `set` if the data support it.
.choose_set <- function(method, std, set = NULL) {
  info <- METHOD_INFO[[method]]
  if (is.null(info)) {
    stop("unknown method '", method, "'.", call. = FALSE)
  }
  avail <- .available_sets(std)
  chosen <- .choose_set_name(method, info, std, set, avail)
  for (v in PREDICTOR_SETS[[chosen]]) {
    na <- which(is.na(std[[v]]))
    if (length(na)) {
      stop("column '", .user_name(v), "' has missing values at ", .row_text(na),
           "; set ", chosen, " needs it complete.", call. = FALSE)
    }
  }
  chosen
}

.user_name <- function(v) {
  map <- c(ITEM_IDS, P9 = PAIN_ID, AGE = "age")
  if (v %in% names(map)) map[[v]] else v
}

.choose_set_name <- function(method, info, std, set, avail) {
  if (!is.null(set)) {
    if (length(set) != 1 || !set %in% info$sets) {
      stop("set '", paste(set, collapse = ","), "' is not a predictor set of ",
           "method '", method, "' (", paste(info$sets, collapse = ", "), ").",
           call. = FALSE)
    }
    if (!set %in% avail) {
      miss <- setdiff(PREDICTOR_SETS[[set]], names(std))
      stop("set ", set, " is not supported by the data; missing columns: ",
           paste(miss, collapse = ", "), ".", call. = FALSE)
    }
    return(set)
  }
  hit <- info$sets[info$sets %in% avail]
  if (length(hit)) return(hit[1])
  if (length(info$sets) == 1) {
    stop("method '", method, "' needs the T-scores, age and sex (set ",
         info$sets, "); missing: ",
         paste(setdiff(PREDICTOR_SETS[[info$sets]], names(std)), collapse = ", "),
         ".", call. = FALSE)
  }
  stop("the data support none of the predictor sets of method '", method,
       "' (", paste(info$sets, collapse = ", "), "); see p2e_columns().",
       call. = FALSE)
}
