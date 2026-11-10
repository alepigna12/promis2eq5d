#!/usr/bin/env Rscript
# The training driver: stages cv, select, bundle, collect (and list).
#
#   Rscript data-raw/train/train.R <stage> --out <dir> --tscore-table <csv>
#       [--smoke] [--commit <sha>] [--dirty 0|1] [unit ids...]
#
# Run on the installed package. Single process; set OPENBLAS_NUM_THREADS=1
# OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 in the shell. The T-score table is a CSV
# with columns domain, raw, tscore; it is read, never written or printed.

args <- commandArgs(trailingOnly = TRUE)
file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
here <- if (length(file_arg)) dirname(normalizePath(sub("^--file=", "", file_arg[1]))) else "data-raw/train"
source(file.path(here, "units.R"))

# ---- arguments ---------------------------------------------------------------
parse_args <- function(a) {
  if (!length(a)) stop("usage: train.R <list|cv|select|bundle|collect> --out <dir> ...", call. = FALSE)
  o <- list(stage = a[1], out = NULL, table = NULL, smoke = FALSE, commit = NA_character_,
            dirty = NA, ids = character())
  i <- 2L
  while (i <= length(a)) {
    k <- a[i]
    if (k == "--smoke") { o$smoke <- TRUE; i <- i + 1L; next }
    if (k %in% c("--out", "--tscore-table", "--commit", "--dirty")) {
      if (i == length(a)) stop(k, " needs a value", call. = FALSE)
      v <- a[i + 1L]
      if (k == "--out") o$out <- v
      if (k == "--tscore-table") o$table <- v
      if (k == "--commit") o$commit <- v
      if (k == "--dirty") {
        if (!v %in% c("0", "1")) stop("--dirty must be 0 or 1, not '", v, "'", call. = FALSE)
        o$dirty <- v == "1"
      }
      i <- i + 2L; next
    }
    if (startsWith(k, "--")) stop("unknown option ", k, call. = FALSE)
    o$ids <- c(o$ids, k); i <- i + 1L
  }
  o
}
opt <- parse_args(args)
SMOKE <- opt$smoke
K_FOLDS <- if (SMOKE) 5L else 10L   # smoke K = 5, not 3: D03 lasso variants run cv.glmnet on K - 1 inner folds, which needs >= 4
Q_DRAWS <- if (SMOKE) 5L else 50L
boot_b <- function(method) {
  r09 <- identical(method, "R09_chain")
  if (SMOKE) { if (r09) 2L else 3L } else { if (r09) 50L else 200L }
}

need_out <- function() if (is.null(opt$out)) stop("--out is required", call. = FALSE)

# ---- helpers -----------------------------------------------------------------
# every write goes to <file>.tmp and is then renamed
write_rds_atomic <- function(x, file, compress = "gzip") {
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  tmp <- paste0(file, ".tmp")
  saveRDS(x, tmp, compress = compress)
  if (!file.rename(tmp, file)) stop("could not rename ", tmp, " to ", file, call. = FALSE)
  invisible(file)
}
write_csv_atomic <- function(x, file) {
  tmp <- paste0(file, ".tmp")
  utils::write.csv(x, tmp, row.names = FALSE)
  if (!file.rename(tmp, file)) stop("could not rename ", tmp, " to ", file, call. = FALSE)
  invisible(file)
}

# Run expr, collecting its warnings (muffled); returns list(value, warnings, seconds)
timed <- function(expr) {
  w <- character()
  t0 <- proc.time()[["elapsed"]]
  v <- withCallingHandlers(expr, warning = function(c) {
    w <<- c(w, conditionMessage(c))
    invokeRestart("muffleWarning")
  })
  list(value = v, warnings = w, seconds = proc.time()[["elapsed"]] - t0)
}

# Variant as a JSON object (scalars only), without a dependency
variant_json <- function(v) {
  one <- function(x) {
    if (is.character(x)) paste0("\"", gsub("\"", "\\\\\"", x), "\"")
    else if (is.logical(x)) tolower(as.character(x))
    else format(x, digits = 17, trim = TRUE)
  }
  paste0("{", paste0("\"", names(v), "\":", vapply(v, one, ""), collapse = ","), "}")
}

# The training data with util_band and (when needed) T-scores as extra
# columns; the stored columns are untouched (the bundle digest covers them
# only). Smoke: a stratified 600-row subsample (country x util_band, in
# proportion), in the original row order.
load_train <- function(need_tscores) {
  suppressPackageStartupMessages(library(promis2eq5d))
  train <- as.data.frame(promis2eq5d:::.training_data())
  train$util_band <- promis2eq5d:::.util_band(as.matrix(train[promis2eq5d:::DIMS]), train$country)
  if (SMOKE) {
    n <- nrow(train)
    stratum <- paste(train$country, train$util_band)
    key <- numeric(n)
    promis2eq5d:::.with_fit_seed(TRAIN_SEED, {
      for (idx in split(seq_len(n), stratum)) {
        len <- length(idx)
        key[idx[sample.int(len)]] <- (seq_len(len) - 0.5) / len   # random rank within the stratum
      }
    })
    train <- train[sort(order(key)[seq_len(600L)]), , drop = FALSE]
    rownames(train) <- NULL
  }
  if (need_tscores) {
    if (is.null(opt$table)) stop("--tscore-table is required for this stage", call. = FALSE)
    tb <- utils::read.csv(opt$table, stringsAsFactors = FALSE)
    tsc <- promis2eq5d:::.tscores_internal(train, tb)
    train[names(tsc)] <- tsc
  }
  train
}

# Folds by country x util_band, from the row order
train_folds <- function(train) {
  promis2eq5d:::make_cv_folds(train, v = K_FOLDS, seed = TRAIN_SEED)$fold
}

# ---- stages ------------------------------------------------------------------
stage_list <- function() {
  what <- opt$ids[1]
  if (is.na(what)) what <- ""
  ids <- switch(what, cv = cv_unit_ids(SMOKE), bundle = bundle_unit_ids(),
                stop("list needs cv or bundle", call. = FALSE))
  cat(ids, sep = "\n")
}

stage_cv <- function() {
  need_out()
  all_ids <- cv_unit_ids(SMOKE)
  ids <- opt$ids
  if (!length(ids)) ids <- all_ids[!file.exists(cv_file(opt$out, all_ids))]
  bad <- setdiff(ids, all_ids)
  if (length(bad)) stop("not a cv unit id: ", paste(bad, collapse = ", "), call. = FALSE)
  train <- load_train(TRUE)
  fold_id <- train_folds(train)
  for (id in ids) {
    u <- parse_unit_id(id)
    m <- promis2eq5d:::.p2e_method(u$method)
    variant <- promis2eq5d:::.p2e_method_variants(m, u$set)[[u$variant_id]]
    r <- timed(promis2eq5d:::.cv_unit(m, train, u$set, u$target, variant, fold_id, TRAIN_SEED))
    res <- c(r$value, list(unit = id, method = u$method, set = u$set, target = u$target,
                           variant_id = u$variant_id, seconds = r$seconds,
                           warnings = r$warnings, smoke = SMOKE))
    write_rds_atomic(res, cv_file(opt$out, id))
    cat(sprintf("%s ok=%s mae=%s %.1fs warnings=%d\n", id, res$ok,
                format(res$mae, digits = 6), r$seconds, length(r$warnings)))
    if (!isTRUE(res$ok)) cat("  reason:", res$reason, "\n")
  }
}

# One row per cell: the selected variant (single-variant cells: variant 1)
selected_table <- function() {
  cells <- train_cells()
  nv <- cell_n_variants(cells, SMOKE)
  ids <- cv_unit_ids(SMOKE)
  miss <- ids[!file.exists(cv_file(opt$out, ids))]
  if (length(miss)) {
    stop("select: ", length(miss), " cv unit(s) missing:\n  ", paste(miss, collapse = "\n  "), call. = FALSE)
  }
  rows <- vector("list", nrow(cells))
  for (i in seq_len(nrow(cells))) {
    m <- promis2eq5d:::.p2e_method(cells$method[i])
    variants <- promis2eq5d:::.p2e_method_variants(m, cells$set[i])
    if (nv[i] == 1L) {
      sel <- data.frame(variant_id = 1L, mae = NA_real_, tie_metric = NA_real_, n_candidates = 1L, n_tied = 1L)
    } else {
      cid <- paste("cv", cells$method[i], cells$set[i], target_label(cells$target[i]), seq_len(nv[i]), sep = ":")
      units <- do.call(rbind, lapply(cid, function(id) {
        r <- readRDS(cv_file(opt$out, id))
        if (!identical(r$smoke, SMOKE)) stop(id, ": the output was made with smoke = ", r$smoke, call. = FALSE)
        data.frame(variant_id = r$variant_id, ok = r$ok, mae = r$mae, tie_metric = r$tie_metric)
      }))
      sel <- promis2eq5d:::.select_variant(units, cell = cell_label(cells$method[i], cells$set[i], cells$target[i]))
      sel <- sel[c("variant_id", "mae", "tie_metric", "n_candidates", "n_tied")]
    }
    rows[[i]] <- data.frame(method = cells$method[i], set = cells$set[i], target = cells$target[i],
                            variant_id = sel$variant_id,
                            variant_json = variant_json(variants[[sel$variant_id]]),
                            mae = sel$mae, tie_metric = sel$tie_metric,
                            n_candidates = sel$n_candidates, n_tied = sel$n_tied,
                            smoke = SMOKE, stringsAsFactors = FALSE)
  }
  do.call(rbind, rows)
}

stage_select <- function() {
  need_out()
  s <- selected_table()
  write_csv_atomic(s, file.path(opt$out, "selected.csv"))
  cat("wrote", file.path(opt$out, "selected.csv"), "with", nrow(s), "cells\n")
}

stage_bundle <- function() {
  need_out()
  if (!SMOKE && (is.na(opt$commit) || !nzchar(opt$commit) || is.na(opt$dirty))) {
    stop("bundle without --smoke needs both --commit <sha> and --dirty 0|1", call. = FALSE)
  }
  options(promis2eq5d.repo_commit = opt$commit, promis2eq5d.repo_dirty = opt$dirty)
  sel <- utils::read.csv(file.path(opt$out, "selected.csv"), stringsAsFactors = FALSE, na.strings = "NA")
  if (!"smoke" %in% names(sel) || !all(sel$smoke == SMOKE)) {
    stop("selected.csv was made with a different smoke setting than this run (smoke = ", SMOKE, ")",
         call. = FALSE)
  }
  all_ids <- bundle_unit_ids()
  ids <- opt$ids
  if (!length(ids)) ids <- all_ids
  bad <- setdiff(ids, all_ids)
  if (length(bad)) stop("not a bundle unit id: ", paste(bad, collapse = ", "), call. = FALSE)
  train <- load_train(TRUE)
  fold_id <- train_folds(train)
  for (id in ids) {
    u <- parse_unit_id(id)
    tmatch <- if (is.na(u$target)) is.na(sel$target) else (!is.na(sel$target) & sel$target == u$target)
    row <- sel[sel$method == u$method & sel$set == u$set & tmatch, ]
    if (nrow(row) != 1L) stop("selected.csv has no unique row for ", id, call. = FALSE)
    m <- promis2eq5d:::.p2e_method(u$method)
    variant <- promis2eq5d:::.p2e_method_variants(m, u$set)[[row$variant_id]]
    file <- promis2eq5d:::.bundle_file(u$method, u$set, u$target)
    B <- boot_b(u$method)
    r <- timed(promis2eq5d:::.fit_bundle(m, train, u$set, u$target, variant, K = K_FOLDS, B = B,
                                         seed = TRAIN_SEED, Q = Q_DRAWS,
                                         fold_id = fold_id))
    write_rds_atomic(r$value, file.path(opt$out, "fits", file), compress = "xz")
    write_rds_atomic(list(unit = id, seconds = r$seconds, warnings = r$warnings, smoke = SMOKE,
                          B = B, K = K_FOLDS, Q = Q_DRAWS),
                     file.path(opt$out, "bundle_logs", file))
    cat(sprintf("%s variant=%d converged=%d attempted=%d failed=%d %.1fs warnings=%d\n", id,
                row$variant_id, length(r$value$boot), r$value$boot_attempted, r$value$boot_failed,
                r$seconds, length(r$warnings)))
  }
}

stage_collect <- function() {
  need_out()
  suppressPackageStartupMessages(library(promis2eq5d))
  cells <- train_cells()
  sp <- file.path(opt$out, "selected.csv")
  s <- if (file.exists(sp)) utils::read.csv(sp, stringsAsFactors = FALSE) else NULL
  rows <- lapply(seq_len(nrow(cells)), function(i) {
    file <- promis2eq5d:::.bundle_file(cells$method[i], cells$set[i], cells$target[i])
    p <- file.path(opt$out, "fits", file)
    vid <- NA_integer_
    if (!is.null(s)) {
      tm <- if (is.na(cells$target[i])) is.na(s$target) else (!is.na(s$target) & s$target == cells$target[i])
      r <- s[s$method == cells$method[i] & s$set == cells$set[i] & tm, ]
      if (nrow(r) == 1L) vid <- r$variant_id
    }
    B <- boot_b(cells$method[i])
    if (!file.exists(p)) {
      return(data.frame(file = file, B = B, converged = NA, short = NA, attempted = NA, failed = NA,
                        variant_id = vid, xz_MB = NA, seconds = NA, warnings = NA))
    }
    b <- readRDS(p)
    lg <- file.path(opt$out, "bundle_logs", file)
    l <- if (file.exists(lg)) readRDS(lg) else list(seconds = NA_real_, warnings = NA)
    data.frame(file = file, B = B, converged = length(b$boot), short = length(b$boot) < B, attempted = b$boot_attempted,
               failed = b$boot_failed, variant_id = vid,
               xz_MB = round(file.size(p) / 1e6, 3), seconds = round(l$seconds, 1),
               warnings = length(l$warnings), stringsAsFactors = FALSE)
  })
  tab <- do.call(rbind, rows)
  print(tab, row.names = FALSE)
  if (any(tab$short %in% TRUE)) {
    cat("
SHORT: converged < B for", sum(tab$short %in% TRUE), "bundle(s):",
        paste(tab$file[tab$short %in% TRUE], collapse = ", "), "
")
  }
  shipped <- cells$method != "R09_chain"
  cat(sprintf("\n%d of %d bundles present; CRAN-shipped bundles (all but R09_chain): %.3f MB (xz)\n",
              sum(!is.na(tab$xz_MB)), nrow(tab), sum(tab$xz_MB[shipped], na.rm = TRUE)))
  cat(sprintf("R09_chain bundles (downloaded): %.3f MB (xz)\n", sum(tab$xz_MB[!shipped], na.rm = TRUE)))
}

switch(opt$stage, list = stage_list(), cv = stage_cv(), select = stage_select(),
       bundle = stage_bundle(), collect = stage_collect(),
       stop("unknown stage '", opt$stage, "'", call. = FALSE))
