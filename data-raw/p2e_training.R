# Builds data/p2e_training.rda from the public Dataverse file (doi:10.7910/DVN/FNYYQR, CC0).
# Usage, from the package directory:
#   Rscript data-raw/p2e_training.R <path to p29_data.sav>
# Ported from iHEA_2027 R/functions/promis.R (orient_items) and R/functions/design.R (load_model_data) at freeze-ihea2027-abstract-v3 (a79e102); changes: only the columns the package uses are kept (no ID, weights, diagnoses or T-scores); P9, items, female and E1-E5 are stored as plain integers.
#
# Orientation: P7_1 to P7_4 are stored with the opposite direction to PROMIS
# scoring and are reversed with 6 - x; every other item is used as stored.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("usage: Rscript data-raw/p2e_training.R <path to p29_data.sav>")
if (!file.exists("DESCRIPTION") || !dir.exists("data-raw")) stop("run from the package directory")
suppressMessages(devtools::load_all(".", quiet = TRUE))

raw <- haven::read_sav(args[[1]])
int <- function(x) as.integer(as.numeric(haven::zap_labels(x)))
need <- c("COUNTRY", "AGE", "GENDER", "P9", ITEMS, DIMS)
if (length(miss <- setdiff(need, names(raw)))) stop("missing columns: ", paste(miss, collapse = ", "))

d <- data.frame(row.names = seq_len(nrow(raw)))
for (v in ITEMS) {
  x <- int(raw[[v]])
  d[[v]] <- if (v %in% paste0("P7_", 1:4)) 6L - x else x
}
d$P9 <- int(raw$P9)
d$AGE <- as.numeric(haven::zap_labels(raw$AGE))
d$female <- as.integer(int(raw$GENDER) == 2L)
d$country <- c("FR", "DE", "UK")[match(int(raw$COUNTRY), 1:3)]
for (v in DIMS) d[[v]] <- int(raw[[v]])
d <- d[c(ITEMS, "P9", "AGE", "female", "country", DIMS)]
for (v in names(d)) attributes(d[[v]]) <- NULL
rownames(d) <- NULL

stopifnot(nrow(d) == 4512L, !anyNA(d),
          all(unlist(d[ITEMS]) %in% 1:5), all(d$P9 %in% 0:10), all(unlist(d[DIMS]) %in% 1:5),
          all(d$female %in% 0:1), all(d$country %in% c("FR", "DE", "UK")),
          all(vapply(d[c(ITEMS, "P9", "female", DIMS)], is.integer, logical(1))),
          is.double(d$AGE), is.character(d$country), !any(grepl("^T_", names(d))))

p2e_training <- d
save(p2e_training, file = "data/p2e_training.rda", compress = "xz", version = 3)
cat("digest:", .data_digest(p2e_training), "\n")
print(table(p2e_training$country))
cat("size (bytes):", file.size("data/p2e_training.rda"), "\n")
