# The cells and units of the training run. Pure functions, no I/O: train.R
# sources this file. A cell is method x predictor set x target; a cv unit is a
# cell x variant (only for cells with more than one variant); a bundle unit is
# a cell; a chunked cell (CHUNKED_METHODS) has a core unit, bootstrap chunk
# units (a range of attempts each) and a merge unit.

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

# The long bundles, fitted as one core unit + bootstrap chunk units + one merge
# unit (stages core, boot, merge) instead of one bundle unit. The others (R09,
# D12) use the bundle stage.
CHUNKED_METHODS <- c("R06_mvop_dim", "D03_clad")
# Default bootstrap attempts per chunk; train.R's --chunk-R06 / --chunk-D03 override them
BOOT_CHUNK_DEFAULT <- c(R06_mvop_dim = 20L, D03_clad = 25L)

chunked_cells <- function() {
  cells <- train_cells()
  cells <- cells[cells$method %in% CHUNKED_METHODS, , drop = FALSE]
  rownames(cells) <- NULL
  cells
}

core_unit_ids <- function() {
  cells <- chunked_cells()
  paste("core", cells$method, cells$set, vapply(cells$target, target_label, ""), sep = ":")
}
merge_unit_ids <- function() {
  cells <- chunked_cells()
  paste("merge", cells$method, cells$set, vapply(cells$target, target_label, ""), sep = ":")
}

# Consecutive ranges "b1-b2" of size `size` over 1..B (the last may be shorter)
boot_ranges <- function(B, size) {
  if (B < 1L) return(character())
  lo <- seq.int(1L, B, by = size)
  paste0(lo, "-", pmin(lo + size - 1L, B))
}

# Bootstrap chunk ids of every chunked cell over 1..B. B: function(method) ->
# B; chunk: chunk size per method (named, as BOOT_CHUNK_DEFAULT)
boot_unit_ids <- function(B, chunk) {
  cells <- chunked_cells()
  unlist(lapply(seq_len(nrow(cells)), function(i) {
    r <- boot_ranges(B(cells$method[i]), chunk[[cells$method[i]]])
    if (!length(r)) return(character())
    paste("boot", cells$method[i], cells$set[i], target_label(cells$target[i]), r, sep = ":")
  }))
}

# "cv:D03_clad:S3:UK:2" -> list(kind, method, set, target (NA for any),
# variant_id, b1, b2); also "bundle|core|merge:<method>:<set>:<target>" and
# "boot:<method>:<set>:<target>:<b1>-<b2>" (b1 and b2 NA for the other kinds)
parse_unit_id <- function(id) {
  p <- strsplit(id, ":", fixed = TRUE)[[1]]
  kind <- p[1]
  ok <- (kind == "cv" && length(p) == 5L) ||
    (kind %in% c("bundle", "core", "merge") && length(p) == 4L) ||
    (kind == "boot" && length(p) == 5L && grepl("^[1-9][0-9]*-[1-9][0-9]*$", p[5]))
  if (!ok) stop("malformed unit id: ", id, call. = FALSE)
  rng <- if (kind == "boot") as.integer(strsplit(p[5], "-", fixed = TRUE)[[1]]) else c(NA_integer_, NA_integer_)
  if (kind == "boot" && rng[1] > rng[2]) stop("malformed unit id (b1 > b2): ", id, call. = FALSE)
  list(kind = kind, method = p[2], set = p[3],
       target = if (p[4] == "any") NA_character_ else p[4],
       variant_id = if (kind == "cv") as.integer(p[5]) else NA_integer_,
       b1 = rng[1], b2 = rng[2])
}

cv_file <- function(out, id) file.path(out, "cv", paste0(gsub(":", "_", id, fixed = TRUE), ".rds"))

# The cell label used in select's messages
cell_label <- function(method, set, target) paste(method, set, target_label(target))
