# Ported from iHEA_2027 R/functions/value_sets.R at freeze-ihea2027-abstract-v3 (a79e102); changes: country codes replaced by p2e_vs value-set objects (eq5d country names via an alias table, or user tables); renamed state_index, score_profiles, vt_min, vt_gap, clamp_utility and vt_additive to .state_index, .score, .vs_min, .vs_gap, .clamp and .vs_additive; dropped the state table, crosswalks and the isotonic fit

# EQ-5D-5L value sets as 3,125-element utility vectors.
#
# State order: element i is the state with levels x = (MO, SC, UA, PD, AD) such
# that i = sum((x - 1) * 5^(4:0)) + 1, so element 1 = 11111, element 2 = 11112
# and element 3125 = 55555 (AD varies fastest).

# Package environment: the alias table (built in .onLoad) and caches
.vs_env <- new.env(parent = emptyenv())
.vs_env$cache <- list()      # p2e_vs objects by code (eq5d value sets only)
.vs_env$additive <- list()   # .vs_additive() results by code (eq5d value sets only)

# eq5d country name, ISO-2 code, ISO-3 code and extra lower-case English
# variants (";"-separated) for each 5L country. Codes of the UK and England
# follow the package ("UK") and ISO 3166-2 ("ENG"). Multi-version countries
# are listed by their base name.
.VS_COUNTRY_TABLE <- data.frame(
  base = c("Australia", "Belgium", "Canada", "China", "Denmark", "Egypt", "England",
           "Ethiopia", "France", "Germany", "Ghana", "HongKong", "Hungary", "India",
           "Indonesia", "Iran", "Ireland", "Italy", "Japan", "Malaysia", "Mexico",
           "Morocco", "Netherlands", "NewZealand", "Nigeria", "Norway", "Peru",
           "Philippines", "Poland", "Portugal", "Romania", "SaudiArabia", "Singapore",
           "Slovenia", "SouthKorea", "Spain", "Sweden", "Taiwan", "Thailand",
           "Trinidad_and_Tobago", "UAE", "Uganda", "UK", "Uruguay", "USA", "Vietnam",
           "WePP"),
  iso2 = c("AU", "BE", "CA", "CN", "DK", "EG", "ENG", "ET", "FR", "DE", "GH", "HK", "HU",
           "IN", "ID", "IR", "IE", "IT", "JP", "MY", "MX", "MA", "NL", "NZ", "NG", "NO",
           "PE", "PH", "PL", "PT", "RO", "SA", "SG", "SI", "KR", "ES", "SE", "TW", "TH",
           "TT", "AE", "UG", "UK", "UY", "US", "VN", NA),
  iso3 = c("AUS", "BEL", "CAN", "CHN", "DNK", "EGY", NA, "ETH", "FRA", "DEU", "GHA", "HKG",
           "HUN", "IND", "IDN", "IRN", "IRL", "ITA", "JPN", "MYS", "MEX", "MAR", "NLD",
           "NZL", "NGA", "NOR", "PER", "PHL", "POL", "PRT", "ROU", "SAU", "SGP", "SVN",
           "KOR", "ESP", "SWE", "TWN", "THA", "TTO", "ARE", "UGA", "GBR", "URY", "USA",
           "VNM", NA),
  variants = c("", "", "", "peoples republic of china", "", "", "", "", "",
               "deutschland", "", "hong kong", "", "", "", "iran;islamic republic of iran",
               "republic of ireland", "", "", "", "", "", "the netherlands;holland",
               "new zealand;nz", "", "", "", "the philippines", "", "", "", "saudi arabia;ksa",
               "", "", "south korea;korea;republic of korea", "", "", "", "", "trinidad and tobago",
               "united arab emirates", "", "united kingdom;gb;great britain;britain", "",
               "united states;united states of america;america", "viet nam", ""),
  stringsAsFactors = FALSE)

# Default eq5d name of each multi-version country (Japan and Peru: the cTTO
# set; Sweden: the most recent set)
.VS_DEFAULTS <- c(Japan = "Japan_cTTO", Peru = "Peru_cTTO", Sweden = "Sweden_2022")

# Alias table (columns alias, eq5d_name) from eq5d::valuesets() and
# .VS_COUNTRY_TABLE, stored in .vs_env$aliases; called from .onLoad()
.build_vs_aliases <- function() {
  vs5 <- eq5d::valuesets(type = "VT", version = "5L")$Country
  base_of <- function(nm) {
    hit <- names(.VS_DEFAULTS)[startsWith(nm, paste0(names(.VS_DEFAULTS), "_"))]
    if (length(hit)) hit[1] else nm
  }
  tab <- .VS_COUNTRY_TABLE
  rows <- list()
  add <- function(alias, nm) {
    alias <- alias[!is.na(alias) & nzchar(alias)]
    if (length(alias)) rows[[length(rows) + 1L]] <<- data.frame(alias = tolower(alias), eq5d_name = nm,
                                                                 stringsAsFactors = FALSE)
  }
  for (nm in vs5) {
    b <- base_of(nm)
    r <- tab[tab$base == b, ]
    if (nrow(r) != 1L) stop("internal error: no country table row for eq5d value set ", nm)
    target <- if (b %in% names(.VS_DEFAULTS)) .VS_DEFAULTS[[b]] else nm
    add(c(nm, gsub("_", " ", nm)), nm)            # full eq5d name, exact
    add(c(b, gsub("_", " ", b), r$iso2, r$iso3, strsplit(r$variants, ";")[[1]]), target)
  }
  al <- do.call(rbind, rows)
  al <- al[!duplicated(al$alias), ]
  rownames(al) <- NULL
  .vs_env$aliases <- al
  invisible(al)
}

# 3,125 x 5 data frame MO...AD of all states in state order. Cached.
.all_states_df <- function() {
  if (is.null(.vs_env$states)) {
    g <- expand.grid(AD = 1:5, PD = 1:5, UA = 1:5, SC = 1:5, MO = 1:5)
    .vs_env$states <- g[, c("MO", "SC", "UA", "PD", "AD")]
    rownames(.vs_env$states) <- NULL
  }
  .vs_env$states
}

# ISO-2-based code of an eq5d name ("UK", "FR", "JP_DCE")
.vs_code <- function(name) {
  b <- names(.VS_DEFAULTS)[startsWith(name, paste0(names(.VS_DEFAULTS), "_"))]
  if (length(b)) {
    paste0(.VS_COUNTRY_TABLE$iso2[.VS_COUNTRY_TABLE$base == b[1]], "_", sub("^[^_]*_", "", name))
  } else {
    .VS_COUNTRY_TABLE$iso2[.VS_COUNTRY_TABLE$base == name]
  }
}

.new_vs <- function(code, name, version, source, u) {
  structure(list(code = code, name = name, version = version, source = source, u = u),
            class = "p2e_vs")
}

# Close alias matches for an unknown country name
.vs_close <- function(x) {
  al <- .vs_env$aliases
  d <- drop(utils::adist(x, al$alias, ignore.case = TRUE))
  ok <- d <= max(2, floor(nchar(x) / 3))
  if (!any(ok)) return(character())
  utils::head(unique(al$eq5d_name[ok][order(d[ok])]), 5)
}

.vs_from_country <- function(country, version) {
  if (!is.character(country) || length(country) != 1L || is.na(country)) {
    stop("country must be a country name or code, a data frame of utilities or a p2e_vs object")
  }
  key <- tolower(trimws(country))
  al <- .vs_env$aliases
  nm <- al$eq5d_name[match(key, al$alias)]
  if (is.na(nm)) {
    cl <- .vs_close(country)
    stop("no EQ-5D-5L value set for '", country, "' in the eq5d package",
         if (length(cl)) paste0(" (close matches: ", paste(cl, collapse = ", "), ")"),
         ". See eq5d::valuesets(type = \"VT\", version = \"5L\").", call. = FALSE)
  }
  b <- names(.VS_DEFAULTS)[startsWith(nm, paste0(names(.VS_DEFAULTS), "_"))]
  suffix <- NA_character_
  code <- NA_character_
  if (length(b)) {
    alts <- unique(al$eq5d_name[startsWith(al$eq5d_name, paste0(b, "_"))])
    sfx <- sub("^[^_]*_", "", alts)
    exact <- key %in% c(tolower(nm), gsub("_", " ", tolower(nm)))   # a full eq5d name was given
    if (!is.null(version)) {
      hit <- which(tolower(sfx) == tolower(as.character(version)))
      if (length(version) != 1L || !length(hit)) {
        stop("no version '", paste(version, collapse = ","), "' for ", b, "; available: ",
             paste(sfx, collapse = ", "), call. = FALSE)
      }
      nm <- alts[hit]; suffix <- sfx[hit]
    } else if (exact) {
      suffix <- sub("^[^_]*_", "", nm)
    } else {
      # the base name or a country code was given: use the documented default and say so
      message("'", country, "' has several EQ-5D-5L value sets (", paste(alts, collapse = ", "),
              "); using ", nm, ". Pass version = (", paste(sfx, collapse = ", "),
              ") to choose another.")
      code <- .VS_COUNTRY_TABLE$iso2[.VS_COUNTRY_TABLE$base == b[1]]
    }
  } else if (!is.null(version)) {
    stop("country '", nm, "' has no alternative versions; version must be NULL", call. = FALSE)
  }
  if (is.na(code)) code <- .vs_code(nm)
  if (is.na(code)) code <- nm                    # WePP has no ISO code
  cached <- .vs_env$cache[[code]]
  if (!is.null(cached)) return(cached)
  u <- eq5d::eq5d(.all_states_df(), version = "5L", type = "VT", country = nm)
  vs <- .new_vs(code, nm, suffix, "eq5d", as.numeric(u))
  .vs_env$cache[[code]] <- vs
  vs
}

.vs_from_table <- function(table) {
  need <- c("MO", "SC", "UA", "PD", "AD", "utility")
  miss <- setdiff(need, names(table))
  if (length(miss)) stop("table must have columns MO, SC, UA, PD, AD and utility; missing: ",
                         paste(miss, collapse = ", "), call. = FALSE)
  if (nrow(table) != 3125L) {
    stop("table must cover all 3,125 EQ-5D-5L states exactly once (it has ", nrow(table), " rows)",
         call. = FALSE)
  }
  idx <- .state_index(as.matrix(table[c("MO", "SC", "UA", "PD", "AD")]))
  if (anyDuplicated(idx)) {
    stop("table must cover all 3,125 EQ-5D-5L states exactly once (duplicated states found)",
         call. = FALSE)
  }
  if (!is.numeric(table$utility) || anyNA(table$utility) || !all(is.finite(table$utility))) {
    stop("table$utility must be numeric and finite", call. = FALSE)
  }
  u <- numeric(3125); u[idx] <- table$utility
  .new_vs("user", "user", NA_character_, "user", u)
}

#' EQ-5D-5L value set from the eq5d package or a user table
#'
#' Creates a value-set object: the 3,125 utilities of the EQ-5D-5L states for
#' one country, taken from the 'eq5d' package, or from a table you supply.
#' Country names, ISO-2 and ISO-3 codes and common English names are accepted.
#'
#' State order: element `i` of `u` belongs to the profile with levels
#' `x = (MO, SC, UA, PD, AD)` with `i = sum((x - 1) * 5^(4:0)) + 1`, so the
#' first element is state 11111 and the last 55555.
#'
#' Japan, Peru and Sweden have several value sets in 'eq5d'. The defaults are
#' `"Japan_cTTO"`, `"Peru_cTTO"` and `"Sweden_2022"`, with a message naming the
#' alternatives; use `version` to choose one.
#'
#' @param country A country name or code (for example `"UK"`, `"gbr"`,
#'   `"United States"`), a data frame with columns `MO`, `SC`, `UA`, `PD`, `AD`
#'   (levels 1 to 5) and `utility` covering all 3,125 states exactly once, or a
#'   value-set object, which is returned unchanged.
#' @param version Optional alternative for countries with several value sets,
#'   as the 'eq5d' name suffix, such as `"DCE"` or `"2020"`.
#'
#' @return An object of class `p2e_vs`: a list with `code` (for example
#'   `"UK"`; `"user"` for a table), `name` (the 'eq5d' country name), `version`
#'   (the alternative suffix, else `NA`), `source` (`"eq5d"` or `"user"`) and
#'   `u`, the 3,125 utilities in state order.
#' @examples
#' uk <- p2e_value_set("UK")
#' length(uk$u)
#' uk$u[1]                       # 11111 is full health
#' p2e_value_set("Japan", version = "DCE")$name
#' @export
p2e_value_set <- function(country, version = NULL) {
  if (inherits(country, "p2e_vs")) return(country)
  if (is.data.frame(country)) return(.vs_from_table(country))
  .vs_from_country(country, version)
}

#' @export
print.p2e_vs <- function(x, ...) {
  cat("<p2e_vs> ", x$name, " (", x$source, "), ", length(x$u), " states, utility range [",
      format(min(x$u), digits = 4), ", ", format(max(x$u), digits = 4), "]\n", sep = "")
  invisible(x)
}

# Row index into the state table of each profile: i = sum((x - 1) * 5^(4:0)) + 1.
# profiles: integer matrix (or data frame) n x 5, or array n x M x 5, with
# levels 1-5 in the order MO, SC, UA, PD, AD. Returns an integer vector n, or
# an n x M matrix for an array. An n x 5 x M array cannot be told apart from
# n x M x 5 when M = 5; the dimensions must be in the documented order.
.state_index <- function(profiles) {
  if (is.data.frame(profiles)) profiles <- as.matrix(profiles)
  d <- dim(profiles)
  if (length(d) == 2 && d[2] == 5) {
    x <- profiles
  } else if (length(d) == 3 && d[3] == 5) {
    x <- matrix(profiles, d[1] * d[2], 5)
  } else {
    stop("profiles must be n x 5 or n x M x 5")
  }
  if (!is.numeric(x) || anyNA(x) || !all(x %in% 1:5)) stop("profile levels must be integers 1-5")
  idx <- as.integer(round(as.vector((x - 1) %*% 5^(4:0)))) + 1L
  if (length(d) == 3) matrix(idx, d[1], d[2]) else idx
}

# A value set, or a list of n value sets (one per row) -> list of the distinct
# value sets and the index of each row's set. Rows sharing a set (same code
# for eq5d sets, identical utilities for user sets) share one entry.
.vs_group <- function(vs, n, what = "profiles") {
  if (inherits(vs, "p2e_vs")) return(list(sets = list(vs), which = rep(1L, n)))
  if (!is.list(vs) || !length(vs) || !all(vapply(vs, inherits, logical(1), "p2e_vs"))) {
    stop("vs must be a p2e_vs value set or a list of them")
  }
  if (length(vs) == 1L) return(list(sets = vs, which = rep(1L, n)))
  if (length(vs) != n) stop("the list of value sets must have length 1 or one per row of ", what)
  sets <- list(); which <- integer(n)
  for (i in seq_len(n)) {
    v <- vs[[i]]; hit <- 0L
    cand <- if (i > 1L) c(which[i - 1L], seq_along(sets)) else seq_along(sets)
    for (k in unique(cand)) {
      s <- sets[[k]]
      if (identical(s$code, v$code) && (s$source == "eq5d" || identical(s$u, v$u))) { hit <- k; break }
    }
    if (hit == 0L) { sets[[length(sets) + 1L]] <- v; hit <- length(sets) }
    which[i] <- hit
  }
  list(sets = sets, which = which)
}

# Utility of each profile on a value set. profiles: as for .state_index().
# vs: one p2e_vs, or a list of n of them (one per row). Returns a numeric
# vector n, or an n x M matrix for an array.
.score <- function(profiles, vs) {
  idx <- .state_index(profiles)
  n <- NROW(idx); m <- NCOL(idx)
  g <- .vs_group(vs, n)
  V <- vapply(g$sets, function(s) s$u, numeric(3125))
  dim(V) <- c(3125L, length(g$sets))
  # Entry r of idx is profile (i, j) with i = (r - 1) %% n + 1, so each row's
  # value set repeats once per draw
  u <- V[cbind(as.vector(idx), rep(g$which, m))]
  if (length(dim(profiles)) == 3) matrix(u, n, m) else u
}

# Lowest utility of a value set (55555 in every set of the pipeline)
.vs_min <- function(vs) min(vs$u)

# Gap below full health: 1 - the best utility below 1
.vs_gap <- function(vs) 1 - max(vs$u[vs$u < 1])

# Limit u to [.vs_min(vs), 1]. u is a vector or an n x M matrix; vs is one
# p2e_vs, or a list with one per element of u or per row of a matrix u.
.clamp <- function(u, vs) {
  if (anyNA(u)) stop("u contains missing values")
  if (inherits(vs, "p2e_vs")) {
    lo <- .vs_min(vs)
  } else {
    if (!length(vs) %in% unique(c(1, NROW(u), length(u)))) {
      stop("vs must have length 1, NROW(u) or length(u)")
    }
    lo <- vapply(vs, .vs_min, numeric(1))
  }
  out <- u                                             # keeps dim, names, dimnames
  out[] <- pmin(pmax(as.vector(u), rep_len(lo, length(u))), 1)
  out
}

# Additive decomposition of a value set (the mean of R05, R06 and R10): the
# least-squares fit of the 3,125 values on an intercept and the dummies of
# levels 2-5 of each dimension, value = intercept - sum_j decrement[j, level_j]
# + remainder. A list of intercept, decrement (5 x 5, dimensions x levels,
# level 1 = 0) and remainder (3,125, state order). An additive set (DE, UK) has
# a least-squares remainder below 1e-12 (floating point), which is set to
# exactly 0, so a method using it can skip the remainder. FR's values are an
# additive model rounded to 3 decimals: |remainder| <= 0.00055 (0.000551).
# Cached by code for eq5d value sets only; user tables are always recomputed.
.vs_additive <- function(vs) {
  use_cache <- identical(vs$source, "eq5d")
  if (use_cache && !is.null(.vs_env$additive[[vs$code]])) return(.vs_env$additive[[vs$code]])
  if (is.null(.vs_env$design)) {
    .vs_env$design <- stats::model.matrix(
      ~ ., data.frame(lapply(.all_states_df(), factor, levels = 1:5)))
  }
  f <- stats::lm.fit(.vs_env$design, vs$u)
  dec <- rbind(0, matrix(-f$coefficients[-1], 4, 5))   # coefficients MO2..MO5, SC2.., per dimension
  r <- unname(f$residuals)
  if (max(abs(r)) < 1e-12) r[] <- 0
  out <- list(intercept = unname(f$coefficients[1]), decrement = unname(t(dec)), remainder = r)
  if (use_cache) .vs_env$additive[[vs$code]] <- out
  out
}
