# The cells and units of the training run. Pure functions, no I/O: train.R
# sources this file. A cell is method x predictor set x target; a cv unit is a
# cell x variant (only for cells with more than one variant); a bundle unit is
# a cell.

TRAIN_SEED <- 20261110L

# Shipped cells, in the order of the brief. target is NA for a response method.
train_cells <- function() {
  grid <- function(method, sets, targets) {
    g <- expand.grid(target = targets, set = sets, stringsAsFactors = FALSE)
    data.frame(method = method, set = g$set, target = g$target, stringsAsFactors = FALSE)
  }
  rbind(
    grid("R06_mvop_dim", c("S1", "S3", "S2", "S4"), NA_character_),
    grid("R09_chain",    c("S3", "S1", "S4", "S2"), NA_character_),
    grid("D03_clad",     c("S3", "S4", "S1", "S2"), c("FR", "DE", "UK")),
    grid("D12_klapproth", "S3K",                    c("FR", "DE", "UK"))
  )
}

# Number of variants of each cell (smoke: R09 configurations 1-2 only). Needs
# the installed package for the variant lists.
cell_n_variants <- function(cells, smoke = FALSE) {
  vapply(seq_len(nrow(cells)), function(i) {
    m <- promis2eq5d:::.p2e_method(cells$method[i])
    n <- length(promis2eq5d:::.p2e_method_variants(m, cells$set[i]))
    if (smoke && identical(cells$method[i], "R09_chain")) n <- min(n, 2L)
    n
  }, integer(1))
}

target_label <- function(target) if (is.na(target)) "any" else target

cv_unit_ids <- function(smoke = FALSE) {
  cells <- train_cells()
  nv <- cell_n_variants(cells, smoke)
  ids <- character()
  for (i in seq_len(nrow(cells))) {
    if (nv[i] > 1L) {
      ids <- c(ids, paste("cv", cells$method[i], cells$set[i], target_label(cells$target[i]),
                          seq_len(nv[i]), sep = ":"))
    }
  }
  ids
}

bundle_unit_ids <- function() {
  cells <- train_cells()
  paste("bundle", cells$method, cells$set, vapply(cells$target, target_label, ""), sep = ":")
}

# "cv:D03_clad:S3:UK:2" -> list(kind, method, set, target (NA for any), variant_id)
parse_unit_id <- function(id) {
  p <- strsplit(id, ":", fixed = TRUE)[[1]]
  kind <- p[1]
  ok <- (kind == "cv" && length(p) == 5L) || (kind == "bundle" && length(p) == 4L)
  if (!ok) stop("malformed unit id: ", id, call. = FALSE)
  list(kind = kind, method = p[2], set = p[3],
       target = if (p[4] == "any") NA_character_ else p[4],
       variant_id = if (kind == "cv") as.integer(p[5]) else NA_integer_)
}

cv_file <- function(out, id) file.path(out, "cv", paste0(gsub(":", "_", id, fixed = TRUE), ".rds"))

# The cell label used in select's messages
cell_label <- function(method, set, target) paste(method, set, target_label(target))
