#!/usr/bin/env Rscript
# The training driver: stages cv, select, bundle, core, boot, merge, collect
# (and list).
#
#   Rscript data-raw/train/train.R <stage> --out <dir> --tscore-table <csv>
#       [--smoke] [--commit <sha>] [--dirty 0|1] [unit ids...]
#   Rscript data-raw/train/train.R list <cv|bundle|core|boot|merge> [--smoke]
#       [--chunk-R06 <n>] [--chunk-D03 <n>]
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
  if (!length(a)) {
    stop("usage: train.R <list|cv|select|bundle|core|boot|merge|collect> --out <dir> ...", call. = FALSE)
  }
  o <- list(stage = a[1], out = NULL, table = NULL, smoke = FALSE, commit = NA_character_,
            dirty = NA, ids = character(), chunk = BOOT_CHUNK_DEFAULT)
  i <- 2L
  while (i <= length(a)) {
    k <- a[i]
    if (k == "--smoke") { o$smoke <- TRUE; i <- i + 1L; next }
    if (k %in% c("--out", "--tscore-table", "--commit", "--dirty", "--chunk-R06", "--chunk-D03")) {
      if (i == length(a)) stop(k, " needs a value", call. = FALSE)
      v <- a[i + 1L]
      if (k == "--out") o$out <- v
      if (k == "--tscore-table") o$table <- v
      if (k == "--commit") o$commit <- v
      if (k == "--dirty") {
        if (!v %in% c("0", "1")) stop("--dirty must be 0 or 1, not '", v, "'", call. = FALSE)
        o$dirty <- v == "1"
      }
      if (k %in% c("--chunk-R06", "--chunk-D03")) {
        if (!grepl("^[1-9][0-9]*$", v)) stop(k, " must be a positive whole number, not '", v, "'", call. = FALSE)
        o$chunk[[if (k == "--chunk-R06") "R06_mvop_dim" else "D03_clad"]] <- as.integer(v)
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
                core = core_unit_ids(), boot = boot_unit_ids(boot_b, opt$chunk),
                merge = merge_unit_ids(),
                stop("list needs cv, bundle, core, boot or merge", call. = FALSE))
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

# Without --smoke, the stages that make bundle pieces need --commit and --dirty;
# they become the bundles' repo_commit and repo_dirty
need_provenance <- function(stage) {
  if (!SMOKE && (is.na(opt$commit) || !nzchar(opt$commit) || is.na(opt$dirty))) {
    stop(stage, " without --smoke needs both --commit <sha> and --dirty 0|1", call. = FALSE)
  }
  options(promis2eq5d.repo_commit = opt$commit, promis2eq5d.repo_dirty = opt$dirty)
}

# selected.csv of this run (an error when it was made in the other smoke mode)
read_selected <- function() {
  sel <- utils::read.csv(file.path(opt$out, "selected.csv"), stringsAsFactors = FALSE, na.strings = "NA")
  if (!"smoke" %in% names(sel) || !all(sel$smoke == SMOKE)) {
    stop("selected.csv was made with a different smoke setting than this run (smoke = ", SMOKE, ")",
         call. = FALSE)
  }
  sel
}

# The cell of a parsed unit id: method object, selected variant, bundle file, B
unit_cell <- function(u, sel, id) {
  tmatch <- if (is.na(u$target)) is.na(sel$target) else (!is.na(sel$target) & sel$target == u$target)
  row <- sel[sel$method == u$method & sel$set == u$set & tmatch, ]
  if (nrow(row) != 1L) stop("selected.csv has no unique row for ", id, call. = FALSE)
  m <- promis2eq5d:::.p2e_method(u$method)
  list(m = m, variant_id = row$variant_id,
       variant = promis2eq5d:::.p2e_method_variants(m, u$set)[[row$variant_id]],
       file = promis2eq5d:::.bundle_file(u$method, u$set, u$target), B = boot_b(u$method))
}

# The ids to run: the given ones (each must be in all_ids), else all_ids
stage_ids <- function(all_ids, what) {
  ids <- opt$ids
  if (!length(ids)) ids <- all_ids
  bad <- setdiff(ids, all_ids)
  if (length(bad)) stop("not a ", what, " unit id: ", paste(bad, collapse = ", "), call. = FALSE)
  ids
}

stage_bundle <- function() {
  need_out()
  need_provenance("bundle")
  sel <- read_selected()
  ids <- stage_ids(bundle_unit_ids(), "bundle")
  train <- load_train(TRUE)
  fold_id <- train_folds(train)
  for (id in ids) {
    u <- parse_unit_id(id)
    cl <- unit_cell(u, sel, id)
    B <- cl$B
    r <- timed(promis2eq5d:::.fit_bundle(cl$m, train, u$set, u$target, cl$variant, K = K_FOLDS, B = B,
                                         seed = TRAIN_SEED, Q = Q_DRAWS,
                                         fold_id = fold_id))
    write_rds_atomic(r$value, file.path(opt$out, "fits", cl$file), compress = "xz")
    write_rds_atomic(list(unit = id, seconds = r$seconds, warnings = r$warnings, smoke = SMOKE,
                          B = B, K = K_FOLDS, Q = Q_DRAWS),
                     file.path(opt$out, "bundle_logs", cl$file))
    cat(sprintf("%s variant=%d converged=%d attempted=%d failed=%d %.1fs warnings=%d\n", id,
                cl$variant_id, length(r$value$boot), r$value$boot_attempted, r$value$boot_failed,
                r$seconds, length(r$warnings)))
  }
}

# ---- chunked bundles: core, boot, merge ----------------------------------------
# A chunked bundle (CHUNKED_METHODS) is one core unit (point and fold fits), boot
# units (bootstrap attempts b1..b2 each) and one merge unit (the replacement
# attempts and the bundle). The merged bundle is identical to what bundle makes
# (bar meta$created). Every piece records smoke, commit, dirty and the variant;
# merge refuses pieces that differ from its own run.
core_file <- function(file) file.path(opt$out, "core", file)
boot_file <- function(file, b1, b2) {
  file.path(opt$out, "boot", paste0(sub("\\.rds$", "", file), "_", b1, "-", b2, ".rds"))
}
# md5 of an object (the uses_oof pool a boot unit used)
obj_md5 <- function(x) {
  f <- tempfile(fileext = ".rds")
  on.exit(unlink(f), add = TRUE)
  saveRDS(x, f, version = 3)
  unname(tools::md5sum(f))
}
# The fields of a piece record that must match this run; returns the mismatches
provenance_mismatch <- function(rec, what, variant_id) {
  now <- list(smoke = SMOKE, commit = opt$commit, dirty = opt$dirty, variant_id = variant_id)
  bad <- names(now)[!vapply(names(now), function(k) identical(rec[[k]], now[[k]]), NA)]
  vapply(bad, function(k) sprintf("%s: %s = %s (this run: %s)", what, k, format(rec[[k]]), format(now[[k]])), "")
}
read_core <- function(file, id) {
  p <- core_file(file)
  if (!file.exists(p)) stop(id, ": the core file ", p, " is missing (run the core unit first)", call. = FALSE)
  readRDS(p)
}

stage_core <- function() {
  need_out()
  need_provenance("core")
  sel <- read_selected()
  ids <- stage_ids(core_unit_ids(), "core")
  train <- load_train(TRUE)
  fold_id <- train_folds(train)
  for (id in ids) {
    u <- parse_unit_id(id)
    cl <- unit_cell(u, sel, id)
    r <- timed(promis2eq5d:::.bundle_core(cl$m, train, u$set, u$target, cl$variant, K = K_FOLDS, B = cl$B,
                                          seed = TRAIN_SEED, Q = Q_DRAWS, fold_id = fold_id))
    write_rds_atomic(list(core = r$value, unit = id, seconds = r$seconds, warnings = r$warnings,
                          smoke = SMOKE, commit = opt$commit, dirty = opt$dirty,
                          variant_id = cl$variant_id, B = cl$B, K = K_FOLDS, Q = Q_DRAWS),
                     core_file(cl$file), compress = "xz")
    cat(sprintf("%s variant=%d %.1fs warnings=%d\n", id, cl$variant_id, r$seconds, length(r$warnings)))
  }
}

stage_boot <- function() {
  need_out()
  need_provenance("boot")
  sel <- read_selected()
  ids <- opt$ids
  if (!length(ids)) ids <- boot_unit_ids(boot_b, opt$chunk)
  cells <- core_unit_ids()
  for (id in ids) {
    u <- parse_unit_id(id)
    if (u$kind != "boot" || !paste("core", u$method, u$set, target_label(u$target), sep = ":") %in% cells) {
      stop("not a boot unit id: ", id, call. = FALSE)
    }
    if (u$b2 > boot_b(u$method)) stop(id, ": attempts beyond B = ", boot_b(u$method), call. = FALSE)
  }
  train <- load_train(TRUE)
  fold_id <- train_folds(train)
  for (id in ids) {
    u <- parse_unit_id(id)
    cl <- unit_cell(u, sel, id)
    pool <- NULL
    vs <- NULL
    if (isTRUE(cl$m$uses_oof)) {
      # the residual pool is the core's
      cr <- read_core(cl$file, id)
      bad <- provenance_mismatch(cr, "core", cl$variant_id)
      if (length(bad)) stop(id, ": the core file does not match this run:\n  ", paste(bad, collapse = "\n  "),
                            call. = FALSE)
      pool <- cr$core$pool
      vs <- cr$core$vs
    }
    r <- timed(promis2eq5d:::.boot_attempts(cl$m, train, u$set, u$target, cl$variant, TRAIN_SEED,
                                            u$b1:u$b2, fold_id, vs = vs, pool = pool))
    write_rds_atomic(list(attempts = r$value, unit = id, b1 = u$b1, b2 = u$b2, seconds = r$seconds,
                          warnings = r$warnings, smoke = SMOKE, commit = opt$commit, dirty = opt$dirty,
                          variant_id = cl$variant_id, pool_md5 = if (!is.null(pool)) obj_md5(pool)),
                     boot_file(cl$file, u$b1, u$b2))
    ok <- sum(!vapply(r$value, function(a) is.null(a$fit), NA))
    cat(sprintf("%s variant=%d converged=%d of %d %.1fs warnings=%d\n", id, cl$variant_id, ok,
                length(r$value), r$seconds, length(r$warnings)))
  }
}

stage_merge <- function() {
  need_out()
  need_provenance("merge")
  sel <- read_selected()
  ids <- stage_ids(merge_unit_ids(), "merge")
  # check every piece before loading the training data
  pieces <- lapply(ids, function(id) {
    u <- parse_unit_id(id)
    cl <- unit_cell(u, sel, id)
    cr <- read_core(cl$file, id)
    stem <- sub("\\.rds$", "", cl$file)
    bf <- list.files(file.path(opt$out, "boot"), full.names = TRUE,
                     pattern = paste0("^", gsub(".", "\\.", stem, fixed = TRUE), "_[0-9]+-[0-9]+\\.rds$"))
    bt <- lapply(bf, readRDS)
    bad <- provenance_mismatch(cr, "core", cl$variant_id)
    for (j in seq_along(bt)) {
      bad <- c(bad, provenance_mismatch(bt[[j]], basename(bf[j]), cl$variant_id))
      if (!identical(basename(bf[j]), basename(boot_file(cl$file, bt[[j]]$b1, bt[[j]]$b2)))) {
        bad <- c(bad, paste0(basename(bf[j]), ": holds attempts ", bt[[j]]$b1, "-", bt[[j]]$b2))
      }
    }
    if (!identical(cr$B, cl$B) || !identical(cr$K, K_FOLDS) || !identical(cr$Q, Q_DRAWS)) {
      bad <- c(bad, sprintf("core: B, K, Q = %d, %d, %d (this run: %d, %d, %d)", cr$B, cr$K, cr$Q,
                            cl$B, K_FOLDS, Q_DRAWS))
    }
    if (isTRUE(cl$m$uses_oof)) {
      pm <- obj_md5(cr$core$pool)
      for (j in seq_along(bt)) {
        if (!identical(bt[[j]]$pool_md5, pm)) bad <- c(bad, paste0(basename(bf[j]), ": made with another core's pool"))
      }
    }
    bs <- unlist(lapply(bt, function(x) x$b1:x$b2))
    dup <- unique(bs[duplicated(bs)])
    miss <- setdiff(seq_len(cl$B), bs)
    out <- setdiff(bs, seq_len(cl$B))
    if (length(dup)) bad <- c(bad, paste("attempts in more than one boot file:", paste(sort(dup), collapse = ", ")))
    if (length(miss)) bad <- c(bad, paste("attempts in no boot file:", paste(miss, collapse = ", ")))
    if (length(out)) bad <- c(bad, paste("attempts beyond B:", paste(sort(out), collapse = ", ")))
    if (length(bad)) stop(id, ": the pieces do not match this run or do not cover 1..", cl$B, ":\n  ",
                          paste(bad, collapse = "\n  "), call. = FALSE)
    bt <- bt[order(vapply(bt, function(x) x$b1, 1L))]
    list(id = id, u = u, cl = cl, core = cr, boot = bt)
  })
  train <- load_train(TRUE)
  for (p in pieces) {
    u <- p$u; cl <- p$cl
    attempts <- unlist(lapply(p$boot, function(x) x$attempts), recursive = FALSE)
    r <- timed(promis2eq5d:::.bundle_assemble(cl$m, train, u$set, u$target, cl$variant, p$core$core,
                                              attempts, cl$B, TRAIN_SEED))
    seconds <- p$core$seconds + sum(vapply(p$boot, function(x) x$seconds, 1)) + r$seconds
    warnings <- c(p$core$warnings, unlist(lapply(p$boot, function(x) x$warnings)), r$warnings)
    write_rds_atomic(r$value, file.path(opt$out, "fits", cl$file), compress = "xz")
    write_rds_atomic(list(unit = p$id, seconds = seconds, warnings = as.character(warnings), smoke = SMOKE,
                          B = cl$B, K = K_FOLDS, Q = Q_DRAWS,
                          parts = c(p$core$unit, vapply(p$boot, function(x) x$unit, ""))),
                     file.path(opt$out, "bundle_logs", cl$file))
    cat(sprintf("%s variant=%d converged=%d attempted=%d failed=%d %.1fs warnings=%d\n", p$id,
                cl$variant_id, length(r$value$boot), r$value$boot_attempted, r$value$boot_failed,
                seconds, length(warnings)))
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
       bundle = stage_bundle(), core = stage_core(), boot = stage_boot(), merge = stage_merge(),
       collect = stage_collect(),
       stop("unknown stage '", opt$stage, "'", call. = FALSE))
