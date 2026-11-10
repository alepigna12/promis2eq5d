# Refit a direct method on the training data rescored under a new value set
# (retrain = TRUE). Response methods rescore under any value set and need no
# refit.

# The training data: the option promis2eq5d.training_data, else the internal
# p2e_training object (shipped by a later plan), else an error.
.training_data <- function() {
  td <- getOption("promis2eq5d.training_data")
  if (!is.null(td)) return(td)
  ns <- asNamespace("promis2eq5d")
  if (exists("p2e_training", envir = ns, inherits = FALSE)) {
    return(get("p2e_training", envir = ns, inherits = FALSE))
  }
  stop("the training data are not available: set options(promis2eq5d.training_data = <data frame>) ",
       "or install a version of promis2eq5d that ships them", call. = FALSE)
}

.retrain_bundle <- function(method, set, vs, boot = NULL, tscore_table = NULL, cache = NULL,
                            interval = FALSE, mode = "mean") {
  m <- .p2e_method(method)
  if (!identical(m$family, "direct")) {
    stop(m$id, " is a response method: it rescores under any value set and needs no refit ",
         "(use retrain = FALSE)", call. = FALSE)
  }
  if (!inherits(vs, "p2e_vs")) stop("vs must be a p2e_vs object", call. = FALSE)
  if (!(is.character(mode) && length(mode) == 1L && mode %in% c("mean", "individual", "cohort"))) {
    stop("mode must be one of \"mean\", \"individual\" or \"cohort\"", call. = FALSE)
  }
  if (!(is.logical(interval) && length(interval) == 1L && !is.na(interval))) {
    stop("interval must be a single TRUE or FALSE", call. = FALSE)
  }
  if (!is.null(cache) && !(is.logical(cache) && length(cache) == 1L && !is.na(cache))) {
    stop("cache must be NULL, TRUE or FALSE", call. = FALSE)
  }

  shipped <- .load_bundle(m$id, set, "UK")
  if (is.null(boot)) boot <- length(shipped$boot)
  if (!(is.numeric(boot) && length(boot) == 1L && !is.na(boot) && boot >= 0)) {
    stop("boot must be a single non-negative number", call. = FALSE)
  }
  B <- if (identical(mode, "mean")) 0L else as.integer(boot)
  K <- max(shipped$fold_id)

  train <- as.data.frame(.training_data())
  miss <- setdiff(.digest_cols(), names(train))
  if (length(miss)) {
    stop("the training data do not match the rows the shipped fits were trained on (missing column(s) ",
         paste(miss, collapse = ", "), "); retrain needs the shipped p2e_training", call. = FALSE)
  }
  digest <- .data_digest(train)
  if (!identical(digest, shipped$meta$train_digest)) {
    stop("the training data do not match the rows the shipped fits were trained on (digest ", digest,
         ", shipped ", shipped$meta$train_digest, "); retrain needs the shipped p2e_training",
         call. = FALSE)
  }
  if (!is_item_set(set)) {
    if (is.null(tscore_table)) {
      stop("retrain on a T-score set needs tscore_table; use an item set (S1 or S2) or retrain = FALSE",
           call. = FALSE)
    }
    tsc <- .tscores_internal(train, tscore_table)
    train[names(tsc)] <- tsc
  }

  tmp <- tempfile(fileext = ".rds")
  on.exit(unlink(tmp), add = TRUE)
  # the key covers the value set, package version, interval, B, the training data
  # actually fitted (after T-scoring) and the T-score table
  hash <- function(x) {
    f <- tempfile(fileext = ".rds"); on.exit(unlink(f), add = TRUE)
    saveRDS(x, f, version = 3)
    unname(tools::md5sum(f))
  }
  saveRDS(list(vs$u, utils::packageVersion("promis2eq5d"), interval, B,
               hash(train), if (is.null(tscore_table)) NULL else hash(tscore_table)),
          tmp, version = 3)
  key <- unname(tools::md5sum(tmp))
  dir <- file.path(.cache_dir(), "retrain")
  file <- file.path(dir, paste0(m$id, "_", set, "_", vs$code, "_", key, ".rds"))

  # reading is not writing: an existing valid file is used unless cache = FALSE
  if (!isFALSE(cache) && file.exists(file)) {
    cached <- tryCatch(readRDS(file), error = function(e) NULL)
    if (inherits(cached, "p2e_bundle") && identical(cached$meta$format, BUNDLE_FORMAT)) {
      message("Using cached refit ", file)
      return(cached)
    }
    if (inherits(cached, "p2e_bundle")) {
      message("The cached refit ", file, " has an old format; it will be refitted")
    } else warning("the cached file '", file, "' is not a readable promis2eq5d fit bundle; ",
            "it will be refitted", call. = FALSE)
  }

  b <- .fit_bundle(m, train, set, vs$code, shipped$variant, K = K, B = B, seed = shipped$meta$seed,
                   progress = TRUE, vs = vs, fold_id = shipped$fold_id, fit_folds = interval)

  save <- isTRUE(cache)
  if (is.null(cache) && interactive()) {
    save <- isTRUE(utils::askYesNo(paste0("Save refitted models to ", dir, " for reuse?")))
  }
  if (save) {
    dir.create(dir, recursive = TRUE, showWarnings = FALSE)
    saveRDS(b, file)
  }
  b
}
