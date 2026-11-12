# Fit bundles: the point fit, K fold fits, B bootstrap fits and out-of-fold
# predictions of one method on one training set (class p2e_bundle).

DATA_DOI <- "10.7910/DVN/FNYYQR"
IHEA_TAG <- "freeze-ihea2027-abstract-v3"

.cache_dir <- function() {
  getOption("promis2eq5d.cache_dir", tools::R_user_dir("promis2eq5d", "cache"))
}

.require_pkgs <- function(pkgs, what) {
  for (p in pkgs) {
    if (!requireNamespace(p, quietly = TRUE)) {
      stop("package '", p, "' is needed to ", what, "; install it with install.packages(\"", p, "\")",
           call. = FALSE)
    }
  }
}

# The bundle format this version reads and writes. Bump it when the stored
# structure changes; .load_bundle() refuses other formats.
BUNDLE_FORMAT <- 2L

# repo_commit and repo_dirty come from options set by the training driver
# (NA when unset): the package does not look for git itself
.bundle_meta <- function(method, seed, K, Q, B, n, train_digest) {
  pk <- unique(c(method$needs, "promis2eq5d", "eq5d"))
  ver <- vapply(pk, function(p) {
    tryCatch(as.character(utils::packageVersion(p)), error = function(e) NA_character_)
  }, character(1))
  list(format = BUNDLE_FORMAT, data_doi = DATA_DOI, seed = seed,
       rng_kind = c("Mersenne-Twister", "Inversion", "Rejection"),
       K = as.integer(K), Q = as.integer(Q), B = as.integer(B), n = as.integer(n),
       train_digest = train_digest,
       repo_commit = getOption("promis2eq5d.repo_commit", NA_character_),
       repo_dirty = getOption("promis2eq5d.repo_dirty", NA),
       ihea_tag = IHEA_TAG, r_version = R.version.string, pkg_versions = ver, created = Sys.time())
}

.digest_cols <- function() c(ITEMS, "P9", "AGE", "female", "country", DIMS)

# The md5 of the training rows a bundle was fitted on: the predictors, covariates,
# country and outcome dimensions as a plain data frame. The data are serialized
# (version 2, xdr) and the 14-byte header (format, writer R version, minimum reader
# version) is dropped; the uncompressed payload is hashed. Unlike the md5 of a
# saveRDS() file, this does not depend on the R version, the OS or zlib.
.data_digest <- function(train) {
  d <- as.data.frame(train)[.digest_cols()]
  rownames(d) <- NULL
  r <- serialize(d, NULL, version = 2, xdr = TRUE)
  if (!identical(r[1:2], charToRaw("X\n")) || !identical(as.integer(r[3:6]), c(0L, 0L, 0L, 2L))) {
    stop("unexpected serialization header", call. = FALSE)
  }
  f <- tempfile(fileext = ".bin")
  on.exit(unlink(f), add = TRUE)
  writeBin(r[-(1:14)], f)
  unname(tools::md5sum(f))
}

# File name of a bundle: <method>_<set>_<target or "any">.rds
.bundle_file <- function(method, set, target) {
  paste0(method, "_", set, "_", if (is.na(target)) "any" else target, ".rds")
}

# A fit did not converge when f$converged is FALSE or f$fit_info$converged is
# not TRUE (the pipeline rule; R06 and R09 report fit_info only). A fit with
# neither field counts as converged.
.fit_converged <- function(f) {
  if (isFALSE(f$converged)) return(FALSE)
  if (is.list(f[["fit_info"]]) && !isTRUE(f[["fit_info"]]$converged)) return(FALSE)
  TRUE
}

# A replicate fails when fit() errors or the fit did not converge
# (.fit_converged). Returns list(fit, reason): the fit and NULL, or NULL and
# the condition message ("did not converge" for a fit that ran, "fit() returned NULL"
# for a fit() that returned NULL).
.try_fit <- function(method, train, set, variant, oof = NULL) {
  err <- NULL
  f <- tryCatch(method$fit(train, set, variant, oof), error = function(e) {
    err <<- conditionMessage(e)
    NULL
  })
  if (is.null(f)) return(list(fit = NULL, reason = if (is.null(err)) "fit() returned NULL" else err))
  if (!.fit_converged(f)) return(list(fit = NULL, reason = "did not converge"))
  list(fit = f, reason = NULL)
}

# withr::with_seed with the RNG kinds pinned, so a seed gives the same stream
# whatever the caller's RNGkind() or RNG version
.with_fit_seed <- function(seed, expr) {
  withr::with_seed(seed, expr, .rng_kind = "Mersenne-Twister", .rng_normal_kind = "Inversion",
                   .rng_sample_kind = "Rejection")
}

# A checksum of a bootstrap index vector, to detect a sample that cannot be rebuilt
.idx_check <- function(idx) sum(as.numeric(idx) * seq_along(idx)) %% 2147483647

# The residual pool of a fit of a uses_oof bundle: the bundle's pool for the
# point fit, else the pool rows of the bootstrap sample of replicate fit$boot_b,
# rebuilt from the seed and the stored countries.
.boot_pool <- function(bundle, fit) {
  if (is.null(fit$boot_b)) return(bundle$pool)
  n <- length(bundle$boot_country)
  idx <- .with_fit_seed(bundle$meta$seed + fit$boot_b, .boot_index(bundle$boot_country, n))
  if (!isTRUE(.idx_check(idx) == fit$idx_check)) {
    stop("the bootstrap sample of replicate ", fit$boot_b, " cannot be rebuilt in this R session",
         call. = FALSE)
  }
  bundle$pool[idx, , drop = FALSE]
}

# Bootstrap row indices: resampled within each country (in radix-sorted order, so the locale does not matter) when
# there is a country column, else simply. Uses the caller's RNG stream.
.boot_index <- function(country, n) {
  if (is.null(country)) return(sample.int(n, n, replace = TRUE))
  unlist(lapply(sort(unique(country), method = "radix"), function(g) {
    w <- which(country == g)
    w[sample.int(length(w), length(w), replace = TRUE)]
  }), use.names = FALSE)
}

# The point and fold fits must fit and converge (as the pipeline's
# .checked_fit): an error from fit() or a fit that did not converge (.fit_converged) stops with the
# method, set and which fit failed
.fit_required <- function(method, set, what, expr) {
  f <- tryCatch(expr, error = function(e) {
    stop(method$id, " (set ", set, "): ", what, " failed: ", class(e)[1], ": ", conditionMessage(e),
         call. = FALSE)
  })
  if (!.fit_converged(f)) {
    stop(method$id, " (set ", set, "): ", what, " did not converge", call. = FALSE)
  }
  f
}

# B: the number of converged bootstrap replicates wanted. Attempts b = 1, 2, ...
# (seed + b) continue until B converge or 2B attempts; fewer than B converged
# warns, none is an error. B = 0 makes no attempts.
# vs: a p2e_vs object used instead of p2e_value_set(target) (direct methods;
#   target becomes vs$code). fold_id: precomputed folds (length nrow(train),
#   values 1..K) instead of make_cv_folds(). fit_folds = FALSE: no fold fits
#   (folds = list(), oof = NULL); a uses_oof method then fits the point model
#   with oof = NULL whatever B is (it computes its own residual pairs, used
#   for the bootstrap pool and for drawing).
# It is .bundle_core() (point and fold fits), .boot_attempts() for attempts
# 1..B and .bundle_assemble() (replacement attempts B + 1, ..., the warning or
# error, the bundle). The training driver runs these pieces in separate tasks
# for the long bundles; the result is identical (bar meta$created).
.fit_bundle <- function(method, train, set, target, variant, K, B, seed, Q = 50,
                        progress = FALSE, vs = NULL, fold_id = NULL, fit_folds = TRUE) {
  core <- .bundle_core(method, train, set, target, variant, K, B, seed, Q = Q, progress = progress,
                       vs = vs, fold_id = fold_id, fit_folds = fit_folds)
  attempts <- .boot_attempts(method, train, set, target, variant, seed, seq_len(B), core$fold_id,
                             vs = core$vs, pool = core$pool, progress = progress, B = B)
  .bundle_assemble(method, train, set, core$target, variant, core, attempts, B, seed, progress = progress)
}

# The value set of a direct method (vs, else p2e_value_set(target)); NULL for a
# response method
.bundle_vs <- function(method, target, vs = NULL) {
  if (!identical(method$family, "direct")) return(NULL)
  if (is.null(vs)) return(p2e_value_set(target))
  if (!inherits(vs, "p2e_vs")) stop("vs must be a p2e_vs object", call. = FALSE)
  vs
}

# The training data every fit of a bundle sees: train with y (direct methods:
# the utility under the value set) and the fold column. Each piece of a bundle
# (.bundle_core, .boot_attempts, .bundle_assemble) recomputes it, so every task
# fits on the same data.
.bundle_train <- function(method, train, target, vs, fold_id) {
  vs <- .bundle_vs(method, target, vs)
  if (!is.null(vs)) train[["y"]] <- .score(as.matrix(train[DIMS]), vs)
  train[["fold"]] <- fold_id
  train
}

# The bootstrap strata: train$country as character (an error when missing)
.boot_strata <- function(train) {
  if (is.null(train[["country"]])) {
    stop("train needs a country column (the bootstrap is stratified by country)", call. = FALSE)
  }
  if (anyNA(train[["country"]])) {
    stop("train$country has missing values; the bootstrap is stratified by country", call. = FALSE)
  }
  as.character(train[["country"]])
}

# A uses_oof replicate or point fit without its residual pairs and fold ids
.drop_oof <- function(f) {
  f$oof_mu <- NULL; f$oof_resid <- NULL; f$foldid <- NULL
  f
}

# Everything of a bundle before the bootstrap: validation (incl. B < 500), the
# value set, the folds, the point and fold fits, oof, the uses_oof residual
# pool, y_profiles and state_freq (arguments as .fit_bundle). Returns a
# p2e_bundle_core: only what the final bundle holds (no training data, no
# T-score columns), plus the settings and data digest .bundle_assemble() checks.
.bundle_core <- function(method, train, set, target, variant, K, B, seed, Q = 50,
                         progress = FALSE, vs = NULL, fold_id = NULL, fit_folds = TRUE) {
  .require_pkgs(method$needs, paste0("fit ", method$id))
  say <- function(...) if (isTRUE(progress)) message(...)
  direct <- identical(method$family, "direct")
  uses_oof <- direct && isTRUE(method$uses_oof)
  n <- nrow(train)
  if (2L * B >= 1000L) {
    stop("B (the number of bootstrap replicates; `boot` in retrain) must be below 500 ",
         "(the attempt seeds must stay below the fold seeds)", call. = FALSE)
  }

  vs <- .bundle_vs(method, target, vs)
  y_profiles <- as.matrix(train[DIMS])
  storage.mode(y_profiles) <- "integer"
  dimnames(y_profiles) <- NULL
  state_freq <- tabulate(.state_index(y_profiles), 3125L)

  if (is.null(fold_id)) {
    fold_id <- make_cv_folds(train, v = K, seed = seed)$fold
    if (anyNA(fold_id) || any(!fold_id %in% seq_len(K)) || any(tabulate(fold_id, K) == 0L)) {
      stop("the computed folds are not valid (values in 1..K, every fold non-empty): ",
           "check that train has the country and util_band columns and enough rows for K = ", K,
           call. = FALSE)
    }
  } else {
    if (length(fold_id) != n || anyNA(fold_id) || any(!fold_id %in% seq_len(K))) {
      stop("fold_id must have length nrow(train) and values in 1..K", call. = FALSE)
    }
    if (any(tabulate(fold_id, K) == 0L)) {
      stop("fold_id must use every fold 1..K (a fold is empty)", call. = FALSE)
    }
  }
  boot_country <- .boot_strata(train)
  train <- .bundle_train(method, train, target, vs, fold_id)

  point <- NULL
  folds <- list()
  oof <- NULL
  oof_pool <- NULL
  point_fit <- function(oof_arg) {
    say("point fit")
    .fit_required(method, set, "point fit",
                  .with_fit_seed(seed, method$fit(train, set, variant, oof_arg)))
  }

  if (!uses_oof) point <- point_fit(NULL)

  if (fit_folds) {
    # oof = FALSE: no residual pairs needed, the fit is used for point predictions only
    fold_oof <- if (uses_oof) FALSE else NULL
    folds <- lapply(seq_len(K), function(k) {
      say("fold ", k, "/", K)
      .fit_required(method, set, paste0("fold ", k, " fit"),
        .with_fit_seed(seed + 1000L + k,
                         method$fit(train[fold_id != k, , drop = FALSE], set, variant, fold_oof)))
    })

    if (direct) {
      oof_mu <- numeric(n)
      for (k in seq_len(K)) {
        te <- fold_id == k
        oof_mu[te] <- method$point(folds[[k]], train[te, , drop = FALSE])
      }
      oof <- list(mu = oof_mu)
      if (uses_oof) oof_pool <- data.frame(oof_mu = oof_mu, oof_resid = train[["y"]] - oof_mu)
    } else {
      prof <- array(NA_integer_, c(n, Q, 5L))
      for (k in seq_len(K)) {
        te <- which(fold_id == k)
        # seeded, so the profiles do not depend on (or advance) the caller's stream
        prof[te, , ] <- .with_fit_seed(seed + 2000L + k,
          method$simulate(folds[[k]], train[te, , drop = FALSE], Q, key = seq_along(te)))
      }
      storage.mode(prof) <- "integer"
      stopifnot(identical(dim(prof), as.integer(c(n, Q, 5L))))
      # the state index of each profile (n x Q), not the n x Q x 5 array
      oof <- list(state = .state_index(prof))
    }
    if (uses_oof) point <- point_fit(oof_pool)
    if (!is.null(method$slim)) folds <- lapply(folds, method$slim)
  } else if (uses_oof) {
    # the point fit computes its own residual pairs whatever B is, so a
    # retrained uses_oof method can always draw (the D03 mean is Monte Carlo
    # over the point fit's draws); the bootstrap pool is those pairs
    point <- point_fit(NULL)
    oof_pool <- data.frame(oof_mu = point$oof_mu, oof_resid = point$oof_resid)
  }
  if (uses_oof) point <- .drop_oof(point)
  train_digest <- .data_digest(train)

  structure(list(
    method = method$id, set = set, target = if (direct) vs$code else NA_character_,
    vs = vs, variant = variant, point = point, folds = folds, fold_id = fold_id,
    pool = oof_pool, boot_country = if (uses_oof) boot_country,
    oof = oof, y_profiles = y_profiles, state_freq = state_freq,
    K = K, Q = Q, B = B, seed = seed, n = n, train_digest = train_digest
  ), class = "p2e_bundle_core")
}

# Fits bootstrap attempts b (an integer vector: any subset of 1..999, any
# order) as .fit_bundle() does: attempt b resamples train within country under
# seed + b (.boot_index) and fits it (.try_fit); a uses_oof method gets
# pool[idx, ] (the core's $pool: required), and its fit loses its residual
# pairs and keeps boot_b and idx_check; slim() is applied. fold_id and vs are
# the core's (fold_id is needed: the fits see train$fold). Needs no core
# otherwise, so attempts of a method without a pool can be fitted before or
# alongside the core. Returns one list(b, fit, reason) per attempt, in the order
# of b (fit NULL and the reason on failure). B and converged only feed the
# progress message "boot attempt b (k of B converged)", k = converged plus the
# attempts converged so far in this call.
.boot_attempts <- function(method, train, set, target, variant, seed, b, fold_id, vs = NULL,
                           pool = NULL, progress = FALSE, B = NULL, converged = 0L) {
  .require_pkgs(method$needs, paste0("fit ", method$id))
  direct <- identical(method$family, "direct")
  uses_oof <- direct && isTRUE(method$uses_oof)
  n <- nrow(train)
  if (!is.numeric(b) || anyNA(b) || any(b != round(b)) || any(b < 1 | b >= 1000) || anyDuplicated(b)) {
    stop("bootstrap attempts b must be distinct whole numbers in 1..999", call. = FALSE)
  }
  b <- as.integer(b)
  if (length(fold_id) != n || anyNA(fold_id)) {
    stop("fold_id must have length nrow(train) and no missing values", call. = FALSE)
  }
  # (no attempts, no pool needed: a uses_oof core without folds and with B = 0 may have an empty pool)
  if (uses_oof && length(b) && (!is.data.frame(pool) || nrow(pool) != n)) {
    stop(method$id, " uses out-of-fold residuals: its bootstrap attempts need the core's pool ",
         "(a data frame with one row per training row)", call. = FALSE)
  }
  boot_country <- .boot_strata(train)
  train <- .bundle_train(method, train, target, vs, fold_id)
  out <- vector("list", length(b))
  for (i in seq_along(b)) {
    bi <- b[i]
    if (isTRUE(progress)) {
      message("boot attempt ", bi, if (!is.null(B)) paste0(" (", converged, " of ", B, " converged)"))
    }
    r <- .with_fit_seed(seed + bi, {
      idx <- .boot_index(boot_country, n)
      fr <- .try_fit(method, train[idx, , drop = FALSE], set, variant,
                     if (uses_oof) pool[idx, , drop = FALSE])
      fr$idx <- idx
      fr
    })
    f <- r$fit
    if (!is.null(f)) {
      if (uses_oof) {
        f <- .drop_oof(f)
        f$boot_b <- bi
        f$idx_check <- .idx_check(r$idx)
      }
      if (!is.null(method$slim)) f <- method$slim(f)
      converged <- converged + 1L
    }
    out[[i]] <- list(b = bi, fit = f, reason = r$reason)
  }
  out
}

# The bundle from its core and the records of attempts 1..B (.boot_attempts;
# each b exactly once, any order): takes the converged fits in b order, makes
# the replacement attempts B + 1, ... while fewer than B converged and fewer
# than 2B were made, then warns (fewer than B converged) or stops (none) as
# .fit_bundle() and returns the p2e_bundle. method, set, target (the core's
# $target: vs$code for a direct method, NA for a response method), variant, B
# and seed and the training data must be the core's. meta is made now (repo_commit and
# repo_dirty from the options, created = now).
.bundle_assemble <- function(method, train, set, target, variant, core, attempts, B, seed,
                             progress = FALSE) {
  if (!inherits(core, "p2e_bundle_core")) stop("core must be a p2e_bundle_core", call. = FALSE)
  if (!identical(core$method, method$id) || !identical(core$set, set)) {
    stop("the core is for ", core$method, " (set ", core$set, "), not ", method$id, " (set ", set, ")",
         call. = FALSE)
  }
  same_target <- if (is.na(core$target)) length(target) == 1L && is.na(target) else
    identical(as.character(target), core$target)
  if (!same_target) stop("the core is for target ", core$target, ", not ", target, call. = FALSE)
  if (!identical(core$variant, variant)) stop("the core was fitted with another variant", call. = FALSE)
  if (!isTRUE(all.equal(core$B, B)) || !isTRUE(all.equal(core$seed, seed))) {
    stop("the core was made with B = ", core$B, " and seed = ", core$seed, ", not B = ", B,
         " and seed = ", seed, call. = FALSE)
  }
  if (!identical(core$n, nrow(train)) || !identical(core$train_digest, .data_digest(train))) {
    stop("the training data differ from those the core was fitted on", call. = FALSE)
  }
  bs <- vapply(attempts, function(a) as.integer(a$b), integer(1))
  dup <- unique(bs[duplicated(bs)])
  if (length(dup)) stop("bootstrap attempt(s) duplicated: ", paste(sort(dup), collapse = ", "), call. = FALSE)
  out <- setdiff(bs, seq_len(B))
  if (length(out)) {
    stop("bootstrap attempt(s) outside 1..", B, ": ", paste(sort(out), collapse = ", "), call. = FALSE)
  }
  miss <- setdiff(seq_len(B), bs)
  if (length(miss)) stop("bootstrap attempt(s) missing: ", paste(miss, collapse = ", "), call. = FALSE)

  boot <- list()
  failures <- data.frame(b = integer(), reason = character())
  add <- function(r) {
    if (is.null(r$fit)) {
      failures <<- rbind(failures, data.frame(b = r$b, reason = r$reason))
    } else {
      boot[[length(boot) + 1L]] <<- r$fit
    }
  }
  for (r in attempts[order(bs)]) add(r)
  b <- as.integer(B)
  while (length(boot) < B && b < 2L * B) {
    b <- b + 1L
    add(.boot_attempts(method, train, set, target, variant, seed, b, core$fold_id, vs = core$vs,
                       pool = core$pool, progress = progress, B = B, converged = length(boot))[[1]])
  }
  if (B > 0L && !length(boot)) {
    stop("all ", b, " bootstrap attempts failed for ", method$id, call. = FALSE)
  }
  if (length(boot) < B) {
    warning("only ", length(boot), " of ", B, " bootstrap replicates converged after ", b,
            " attempts for ", method$id, call. = FALSE)
  }

  structure(list(
    method = core$method, set = core$set, target = core$target,
    vs = core$vs, variant = core$variant, point = core$point, folds = core$folds, fold_id = core$fold_id,
    boot = boot, boot_attempted = b, boot_failed = b - length(boot), boot_failures = failures,
    pool = core$pool, boot_country = core$boot_country,
    oof = core$oof, y_profiles = core$y_profiles, state_freq = core$state_freq,
    meta = .bundle_meta(method, core$seed, core$K, core$Q, core$B, core$n, core$train_digest)
  ), class = "p2e_bundle")
}

# Looks for <method>_<set>_<target or "any">.rds in fit_dir, the installed
# fits and the cache. Downloads only after an interactive yes (R09_chain).
.load_bundle <- function(method, set, target) {
  file <- .bundle_file(method, set, target)
  dirs <- c(getOption("promis2eq5d.fit_dir"),
            system.file("fits", package = "promis2eq5d"),
            file.path(.cache_dir(), "fits"))
  dirs <- dirs[nzchar(dirs)]
  find <- function() {
    for (d in dirs) {
      p <- file.path(d, file)
      if (file.exists(p)) return(p)
    }
    NULL
  }
  p <- find()
  if (is.null(p) && identical(method, "R09_chain")) {
    row <- .manifest_rows(method, set)
    if (nrow(row) == 1L && .is_interactive() &&
        isTRUE(.ask_yes_no(sprintf(paste0("Fit bundle '%s' is not installed. Download it now ",
                                          "(about %.0f MB) with p2e_download()?"),
                                   file, as.numeric(row$bytes) / 1e6)))) {
      p2e_download(method, set, .consented = TRUE)
      p <- find()
    }
  }
  if (is.null(p)) {
    hint <- if (identical(method, "R09_chain")) {
      paste0("; run p2e_download('R09_chain', set = '", set, "')")
    } else ""
    lead <- if (identical(method, "R09_chain") && !nrow(.manifest_rows(method, set))) {
      paste0("no download is available yet for ", method, " ", set, "; ")
    } else ""
    stop(lead, "fit bundle '", file, "' not found (searched: ", paste(dirs, collapse = ", "), ")",
         hint, call. = FALSE)
  }
  b <- readRDS(p)
  if (!inherits(b, "p2e_bundle")) stop("'", p, "' is not a promis2eq5d fit bundle")
  .check_bundle_format(b, file)
  b
}
# Stops unless the bundle has the format this version reads
.check_bundle_format <- function(b, file) {
  if (!identical(b$meta$format, BUNDLE_FORMAT)) {
    fmt <- if (is.null(b$meta$format)) "none" else format(b$meta$format)
    stop("fit bundle '", file, "' has format ", fmt, "; this version of promis2eq5d needs format ",
         BUNDLE_FORMAT, ": reinstall promis2eq5d, run p2e_download() again, or replace the file in the folder set by options(promis2eq5d.fit_dir)",
         call. = FALSE)
  }
  invisible(b)
}

