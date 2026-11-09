# Ported from iHEA_2027 R/functions/promis.R at freeze-ihea2027-abstract-v3 (a79e102); changes: PROMIS_DOMAINS only; no orientation step

# Item -> domain. Sleep disturbance is P5 plus P6_1-P6_3. The items are used
# as supplied: they must already be in PROMIS scoring orientation.
PROMIS_DOMAINS <- list(
  PF  = paste0("P1_", 1:4),
  ANX = paste0("P2_", 1:4),
  DEP = paste0("P3_", 1:4),
  FAT = paste0("P4_", 1:4),
  SLP = c("P5", paste0("P6_", 1:3)),
  SR  = paste0("P7_", 1:4),
  PI  = paste0("P8_", 1:4)
)

#' PROMIS T-scores from a conversion table you supply
#'
#' Computes the seven PROMIS-29 domain T-scores from the 28 item responses.
#' Each domain raw score is the sum of its items (4 to 20), mapped to a
#' T-score through `table`. Items must already be in PROMIS scoring orientation,
#' where a higher value means the PROMIS direction of its domain, as in the
#' HealthMeasures scoring manual; the function does not reverse any item.
#'
#' No conversion table is shipped with this package, because PROMIS scoring
#' tables are copyrighted and HealthMeasures requires written permission to
#' redistribute them. Get T-scores from the HealthMeasures Scoring Service, the
#' HealthMeasures API or REDCap, or supply your own licensed table here.
#'
#' @param data A data frame with one row per respondent and the 28 item columns
#'   named by their HealthMeasures IDs (see [p2e_columns()]). Items are
#'   integers 1-5. A missing item gives `NA` for its domain.
#' @param table A data frame with columns `domain` (`"PF"`, `"ANX"`, `"DEP"`,
#'   `"FAT"`, `"SLP"`, `"SR"`, `"PI"`), `raw` (the raw sum, 4-20) and `tscore`
#'   (the T-score for that domain and raw sum).
#' @return A data frame with the seven columns `T_PF`, `T_ANX`, `T_DEP`,
#'   `T_FAT`, `T_SLP`, `T_SR` and `T_PI`, one row per row of `data` in the same
#'   order. A domain with a missing item is `NA`.
#' @export
#' @examples
#' # A made-up table, not PROMIS values: T-score = 50 + raw sum
#' domains <- c("PF", "ANX", "DEP", "FAT", "SLP", "SR", "PI")
#' tb <- expand.grid(domain = domains, raw = 4:20, stringsAsFactors = FALSE)
#' tb$tscore <- 50 + tb$raw
#' cols <- subset(p2e_columns(), type == "item")$column
#' items <- as.data.frame(setNames(as.list(rep(3L, length(cols))), cols))
#' p2e_tscores(items, tb)
p2e_tscores <- function(data, table) {
  if (!is.data.frame(data)) stop("`data` must be a data frame.", call. = FALSE)
  id <- unname(ITEM_IDS)
  miss <- setdiff(id, names(data))
  if (length(miss)) {
    stop("missing item columns: ", paste(miss, collapse = ", "), ".", call. = FALSE)
  }
  std <- .standardise_input(data[id])
  .tscores_internal(std, table)
}

# Compute the seven T-score columns from internal item names (P1_1 ...).
# `std` is a data frame with the 28 ITEMS columns. Returns a data frame with
# the TSCORES columns, one row per row of `std`.
.tscores_internal <- function(std, table) {
  need <- c("domain", "raw", "tscore")
  miss <- setdiff(need, names(table))
  if (length(miss)) {
    stop("the T-score table needs the columns ", paste(need, collapse = ", "),
         "; missing: ", paste(miss, collapse = ", "), ".", call. = FALSE)
  }
  table <- as.data.frame(table, stringsAsFactors = FALSE)
  table$domain <- as.character(table$domain)
  if (!is.numeric(table$raw) || !is.numeric(table$tscore)) {
    stop("the columns 'raw' and 'tscore' of the T-score table must be numeric.",
         call. = FALSE)
  }
  out <- data.frame(row.names = seq_len(nrow(std)))
  for (dm in names(PROMIS_DOMAINS)) {
    raw <- unname(rowSums(std[PROMIS_DOMAINS[[dm]]]))  # NA if any item is NA
    tbl <- table[table$domain == dm, , drop = FALSE]
    bad <- !is.na(raw) & is.na(match(raw, tbl$raw))
    if (any(bad)) {
      stop("no table row for domain ", dm, " at raw score(s) ",
           paste(sort(unique(raw[bad])), collapse = ", "), ".", call. = FALSE)
    }
    out[[paste0("T_", dm)]] <- tbl$tscore[match(raw, tbl$raw)]
  }
  out
}
