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

.git_commit <- function() {
  tryCatch({
    src <- find.package("promis2eq5d")
    out <- suppressWarnings(system2("git", c("-C", shQuote(src), "rev-parse", "HEAD"),
                                    stdout = TRUE, stderr = FALSE))
    if (length(out) == 1L && grepl("^[0-9a-f]{40}$", out)) out else NA_character_
  }, error = function(e) NA_character_, warning = function(w) NA_character_)
}

.bundle_meta <- function(method, seed) {
  pk <- unique(c(method$needs, "promis2eq5d", "eq5d"))
  ver <- vapply(pk, function(p) {
    tryCatch(as.character(utils::packageVersion(p)), error = function(e) NA_character_)
  }, character(1))
  list(data_doi = DATA_DOI, seed = seed, repo_commit = .git_commit(), ihea_tag = IHEA_TAG,
       r_version = R.version.string, pkg_versions = ver, created = Sys.time())
}

# A replicate fails when fit() errors or reports isFALSE(converged)
.try_fit <- function(method, train, set, variant, oof = NULL) {
  f <- tryCatch(method$fit(train, set, variant, oof), error = function(e) NULL)
  if (is.null(f) || isFALSE(f$converged)) NULL else f
}

# Bootstrap row indices: resampled within each country (in sorted order) when
# there is a country column, else simply. Uses the caller's RNG stream.
.boot_index <- function(country, n) {
  if (is.null(country)) return(sample.int(n, n, replace = TRUE))
  unlist(lapply(sort(unique(country)), function(g) {
    w <- which(country == g)
    w[sample.int(length(w), length(w), replace = TRUE)]
  }), use.names = FALSE)
}

# The point and fold fits must fit and converge (as the pipeline's
# .checked_fit): an error from fit() or converged = FALSE stops with the
# method, set and which fit failed
.fit_required <- function(method, set, what, expr) {
  f <- tryCatch(expr, error = function(e) {
    stop(method$id, " (set ", set, "): ", what, " failed: ", conditionMessage(e), call. = FALSE)
  })
  if (isFALSE(f$converged)) {
    stop(method$id, " (set ", set, "): ", what, " did not converge", call. = FALSE)
  }
  f
}

.fit_bundle <- function(method, train, set, target, variant, K, B, seed, Q = 50,
                        progress = FALSE) {
  .require_pkgs(method$needs, paste0("fit ", method$id))
  say <- function(...) if (isTRUE(progress)) message(...)
  direct <- identical(method$family, "direct")
  uses_oof <- direct && isTRUE(method$uses_oof)
  n <- nrow(train)

  vs <- NULL
  if (direct) {
    vs <- p2e_value_set(target)
    train$y <- .score(as.matrix(train[DIMS]), vs)
  }
  y_profiles <- as.matrix(train[DIMS])
  storage.mode(y_profiles) <- "integer"
  dimnames(y_profiles) <- NULL
  state_freq <- tabulate(.state_index(y_profiles), 3125L)

  fold_id <- make_cv_folds(train, v = K, seed = seed)$fold
  train$fold <- fold_id

  point <- NULL
  if (!uses_oof) {
    say("fitting the point model")
    point <- .fit_required(method, set, "point fit",
                           withr::with_seed(seed, method$fit(train, set, variant, NULL)))
  }

  say("fitting ", K, " fold models")
  # oof = FALSE: no residual pairs needed, the fit is used for point predictions only
  fold_oof <- if (uses_oof) FALSE else NULL
  folds <- lapply(seq_len(K), function(k) {
    .fit_required(method, set, paste0("fold ", k, " fit"),
      withr::with_seed(seed + 1000L + k,
                       method$fit(train[fold_id != k, , drop = FALSE], set, variant, fold_oof)))
  })

  say("computing out-of-fold predictions")
  oof_pool <- NULL
  if (direct) {
    oof_mu <- numeric(n)
    for (k in seq_len(K)) {
      te <- fold_id == k
      oof_mu[te] <- method$point(folds[[k]], train[te, , drop = FALSE])
    }
    oof <- list(mu = oof_mu)
    if (uses_oof) oof_pool <- data.frame(oof_mu = oof_mu, oof_resid = train$y - oof_mu)
  } else {
    prof <- array(NA_integer_, c(n, Q, 5L))
    for (k in seq_len(K)) {
      te <- which(fold_id == k)
      prof[te, , ] <- method$simulate(folds[[k]], train[te, , drop = FALSE], Q,
                                      key = seq_along(te))
    }
    oof <- list(profiles = prof)
  }

  if (uses_oof) {
    say("fitting the point model")
    point <- .fit_required(method, set, "point fit",
                           withr::with_seed(seed, method$fit(train, set, variant, oof_pool)))
  }

  boot <- list()
  failed <- 0L
  for (b in seq_len(B)) {
    say("bootstrap replicate ", b, " of ", B)
    f <- withr::with_seed(seed + b, {
      idx <- .boot_index(train$country, n)
      .try_fit(method, train[idx, , drop = FALSE], set, variant,
               if (uses_oof) oof_pool[idx, , drop = FALSE])
    })
    if (is.null(f)) failed <- failed + 1L else boot[[length(boot) + 1L]] <- f
  }
  if (B > 0L && !length(boot)) stop("all ", B, " bootstrap replicates failed for ", method$id)

  structure(list(
    method = method$id, set = set, target = if (direct) vs$code else NA_character_,
    vs = vs, variant = variant, point = point, folds = folds, fold_id = fold_id,
    boot = boot, boot_attempted = as.integer(B), boot_failed = as.integer(failed),
    oof = oof, y_profiles = y_profiles, state_freq = state_freq,
    meta = .bundle_meta(method, seed)
  ), class = "p2e_bundle")
}

# Looks for <method>_<set>_<target or "any">.rds in fit_dir, the installed
# fits and the cache. Never downloads.
.load_bundle <- function(method, set, target) {
  file <- paste0(method, "_", set, "_", if (is.na(target)) "any" else target, ".rds")
  dirs <- c(getOption("promis2eq5d.fit_dir"),
            system.file("fits", package = "promis2eq5d"),
            file.path(.cache_dir(), "fits"))
  dirs <- dirs[nzchar(dirs)]
  for (d in dirs) {
    p <- file.path(d, file)
    if (file.exists(p)) {
      b <- readRDS(p)
      if (!inherits(b, "p2e_bundle")) stop("'", p, "' is not a promis2eq5d fit bundle")
      return(b)
    }
  }
  hint <- if (identical(method, "R09_chain")) {
    paste0("; run p2e_download('R09_chain', set = '", set, "')")
  } else ""
  stop("fit bundle '", file, "' not found (searched: ", paste(dirs, collapse = ", "), ")", hint,
       call. = FALSE)
}
