test_that("a bundle has K fold fits, B boot fits and aligned oof", {
  b <- .fit_bundle(toy_method(), train = make_synthetic_data(200, 1), set = "S4", target = "UK",
                   variant = list(), K = 3, B = 4, seed = 1)
  expect_s3_class(b, "p2e_bundle"); expect_length(b$folds, 3); expect_length(b$boot, 4)
  expect_equal(length(b$oof$mu), 200); expect_equal(sort(unique(b$fold_id)), 1:3)
  expect_equal(sum(b$state_freq), 200); expect_equal(dim(b$y_profiles), c(200L, 5L))
  expect_s3_class(b$vs, "p2e_vs"); expect_equal(b$target, "UK")
})
test_that("bundles carry format 2 and the full meta", {
  withr::local_options(promis2eq5d.repo_commit = "abc", promis2eq5d.repo_dirty = FALSE)
  tr <- make_synthetic_data(200, 1)
  b <- .fit_bundle(toy_method(), tr, "S4", "UK", list(), 3, 2, 1)
  expect_named(b$meta, c("format","data_doi","seed","rng_kind","K","Q","B","n","train_digest",
                         "repo_commit","repo_dirty","ihea_tag","r_version","pkg_versions","created"))
  expect_identical(b$meta$format, 2L); expect_equal(b$meta$n, 200L); expect_equal(b$meta$K, 3L)
  expect_equal(b$meta$Q, 50L); expect_equal(b$meta$B, 2L)
  expect_identical(b$meta$rng_kind, c("Mersenne-Twister", "Inversion", "Rejection"))
  expect_identical(b$meta$train_digest, .data_digest(tr)); expect_identical(b$meta$repo_commit, "abc")
  expect_identical(b$meta$repo_dirty, FALSE)
})
test_that("the data digest ignores row names and other columns, and tracks the data", {
  tr <- make_synthetic_data(60, 1)
  d <- .data_digest(tr)
  expect_match(d, "^[0-9a-f]{32}$")
  t2 <- tr; rownames(t2) <- paste0("r", seq_len(60)); t2$WEIGHT <- 0
  expect_identical(.data_digest(t2), d)
  t3 <- tr; t3$AGE[1] <- t3$AGE[1] + 1L
  expect_false(identical(.data_digest(t3), d))
})
digest_frame <- function() {
  d <- as.data.frame(matrix(rep(c(1L, 2L, 3L), length.out = 3 * length(ITEMS)), 3, length(ITEMS)))
  names(d) <- ITEMS
  d$P9 <- c(0L, 5L, 10L); d$AGE <- c(30.5, 45, 62.25); d$female <- c(0L, 1L, 1L)
  d$country <- c("FR", "DE", "UK")
  for (v in DIMS) d[[v]] <- c(1L, 3L, 5L)
  d
}
test_that("the data digest has a pinned value", {
  expect_identical(.data_digest(digest_frame()), "53ffbf037c7b77f6535e60704a2620d9")
})
test_that("the data digest is the md5 of the serialization payload, without its header", {
  d <- digest_frame()[.digest_cols()]
  f <- tempfile(); on.exit(unlink(f))
  writeBin(serialize(d, NULL, version = 2, xdr = TRUE)[-(1:14)], f)
  expect_identical(.data_digest(d), unname(tools::md5sum(f)))
})
test_that("an old-format bundle is refused with what to do", {
  d <- withr::local_tempdir(); withr::local_options(promis2eq5d.fit_dir = d)
  b <- .fit_bundle(toy_method(), make_synthetic_data(100, 1), "S4", "UK", list(), 3, 0, 1)
  b$meta$format <- NULL; saveRDS(b, file.path(d, .bundle_file("D03_clad", "S4", "UK")))
  expect_error(.load_bundle("D03_clad", "S4", "UK"), "format none.*needs format 2.*reinstall")
})
test_that(".bundle_file names response bundles 'any'", {
  expect_identical(.bundle_file("R09_chain", "S3", NA), "R09_chain_S3_any.rds")
  expect_identical(.bundle_file("D03_clad", "S1", "FR"), "D03_clad_S1_FR.rds")
})
test_that("slim() is applied to fold and bootstrap fits, never to the point fit", {
  m <- toy_method(); m$slim <- function(fit) { fit$slimmed <- TRUE; fit }
  b <- .fit_bundle(m, make_synthetic_data(200, 1), "S4", "UK", list(), 3, 2, 1)
  expect_null(b$point$slimmed)
  for (f in c(b$folds, b$boot)) expect_true(f$slimmed)
})
test_that("a fit() that returns NULL is a replaced replicate with a reason", {
  # calls 1-4 are the point and fold fits; attempt b is call 4 + b
  m <- toy_method(); m$fit <- local({ f <- m$fit; n <- 0L
    function(...) { n <<- n + 1L; if (n == 6L) return(NULL); f(...) } })
  b <- .fit_bundle(m, make_synthetic_data(200, 1), "S4", "UK", list(), 3, 3, 1)
  expect_length(b$boot, 3); expect_equal(b$boot_attempted, 4L)
  expect_equal(b$boot_failures$b, 2L); expect_identical(b$boot_failures$reason, "fit() returned NULL")
})
test_that(".boot_index does not depend on the collation locale", {
  country <- rep(c("b", "B", "a", "A"), length.out = 60)
  # strata in radix (C) order: A, B, a, b, whatever the locale or the test context
  hand <- withr::with_seed(3, unlist(lapply(c("A", "B", "a", "b"), function(g) {
    w <- which(country == g); w[sample.int(length(w), length(w), replace = TRUE)]
  })))
  expect_identical(withr::with_collate("C", withr::with_seed(3, .boot_index(country, 60))), hand)
  loc <- if (.Platform$OS.type == "windows") "English_United States.1252" else "en_US.UTF-8"
  other <- tryCatch(withr::with_collate(loc, sort(c("a", "B"))), error = function(e) NULL, warning = function(w) NULL)
  skip_if(is.null(other) || identical(other, c("B", "a")), "no locale that collates differently from C")
  expect_identical(withr::with_collate(loc, withr::with_seed(3, .boot_index(country, 60))), hand)
})
test_that("progress shows the attempt and how many have converged", {
  m <- toy_method(); m$fit <- local({ f <- m$fit; n <- 0L
    function(...) { n <<- n + 1L; if (n == 5L) stop("flaky"); f(...) } })
  msgs <- testthat::capture_messages(
    .fit_bundle(m, make_synthetic_data(200, 1), "S4", "UK", list(), 3, 2, 1, progress = TRUE))
  expect_equal(msgs[-(1:4)], c("boot attempt 1 (0 of 2 converged)\n", "boot attempt 2 (0 of 2 converged)\n",
                               "boot attempt 3 (1 of 2 converged)\n"))
})
test_that("failed replicates are replaced until B converge", {
  b <- .fit_bundle(toy_method(fail_on = 2), make_synthetic_data(200, 1), "S4", "UK", list(), 3, 4, 1)
  expect_length(b$boot, 4); expect_equal(b$boot_attempted, 5L); expect_equal(b$boot_failed, 1L)
  expect_equal(b$boot_failures$b, 2L); expect_match(b$boot_failures$reason, "failed on purpose")
})
test_that("the attempt cap keeps the converged replicates and warns", {
  # calls 1-4 are the point and the K = 3 fold fits; attempt b is call 4 + b.
  # Only attempts 1, 3 and 5 converge, so the cap of 2B = 8 attempts is reached
  m <- toy_method(); m$fit <- local({ f <- m$fit; n <- 0L
    function(...) { n <<- n + 1L; if (n > 4L && !(n - 4L) %in% c(1L, 3L, 5L)) stop("flaky"); f(...) } })
  expect_warning(b <- .fit_bundle(m, make_synthetic_data(200, 1), "S4", "UK", list(), 3, 4, 1),
                 "3 of 4 .*8 attempts")
  expect_length(b$boot, 3L); expect_equal(b$boot_attempted, 8L)
})
test_that("B = 0 makes no attempts and no warning", {
  expect_no_warning(b <- .fit_bundle(toy_method(), make_synthetic_data(200, 1), "S4", "UK", list(), 3, 0, 1))
  expect_length(b$boot, 0); expect_equal(b$boot_attempted, 0L); expect_equal(b$boot_failed, 0L)
})
test_that("when no replicate converges the fit stops", {
  m <- toy_method(); m$fit <- local({ f <- m$fit; n <- 0L
    function(...) { n <<- n + 1L; if (n > 4L) stop("flaky"); f(...) } })
  expect_error(.fit_bundle(m, make_synthetic_data(200, 1), "S4", "UK", list(), 3, 2, 1), "failed for toy")
})
test_that("B of 500 or more is refused (attempt seeds would reach the fold seeds)", {
  expect_error(.fit_bundle(toy_method(), make_synthetic_data(50, 1), "S4", "UK", list(), 3, 500, 1),
               "B \\(the number of bootstrap replicates; `boot` in retrain\\) must be below 500")
})
test_that("uses_oof bundles store one pool and rebuild each replicate's pool", {
  train <- make_synthetic_data(200, 1)
  b <- .fit_bundle(toy_method(uses_oof = TRUE), train, "S4", "UK", list(), 3, 3, 7)
  expect_named(b$pool, c("oof_mu", "oof_resid")); expect_equal(b$pool$oof_mu, b$oof$mu)
  for (f in b$boot) {
    expect_null(f$oof_mu); expect_true(is.integer(f$boot_b))
    expect_equal(.boot_pool(b, f), f$oof_given, ignore_attr = TRUE)   # the toy keeps what it was given
  }
  expect_equal(.boot_pool(b, b$point), b$pool)
})
test_that("pools are rebuilt under the user's RNG kind", {
  b <- .fit_bundle(toy_method(uses_oof = TRUE), make_synthetic_data(200, 1), "S4", "UK", list(), 3, 2, 7)
  ref <- .boot_pool(b, b$boot[[2]])
  withr::with_rng_version("3.5.0", expect_equal(.boot_pool(b, b$boot[[2]]), ref))
  withr::with_seed(1, .rng_kind = "L'Ecuyer-CMRG", expect_equal(.boot_pool(b, b$boot[[2]]), ref))
})
test_that("a pool that cannot be rebuilt is an error", {
  b <- .fit_bundle(toy_method(uses_oof = TRUE), make_synthetic_data(200, 1), "S4", "UK", list(), 3, 2, 7)
  b$boot[[1]]$idx_check <- b$boot[[1]]$idx_check + 1
  expect_error(.boot_pool(b, b$boot[[1]]), "cannot be rebuilt")
})
test_that("missing bundles give an actionable error", {
  withr::local_options(promis2eq5d.fit_dir = tempfile(), promis2eq5d.cache_dir = tempfile())
  expect_error(.load_bundle("R09_chain", "S3", NA), "p2e_download('R09_chain', set = 'S3')", fixed = TRUE)
})
test_that("a missing Suggests package is an install hint", {
  m <- toy_method(); m$needs <- "notapkg"
  expect_error(.fit_bundle(m, make_synthetic_data(50, 1), "S4", "UK", list(), 3, 2, 1), "notapkg.*install")
})

test_that("oof and bootstrap use the documented folds and seeds", {
  train <- make_synthetic_data(200, 1)
  b <- .fit_bundle(toy_method(), train, "S4", "UK", list(), 3, 2, 7)
  expect_equal(b$fold_id, make_cv_folds(train, v = 3, seed = 7)$fold)
  # oof$mu is the prediction of the fold fit that left the row out
  k <- 2L; te <- b$fold_id == k
  pred <- as.numeric(cbind(1, as.matrix(train[te, PREDICTOR_SETS$S4])) %*% b$folds[[k]]$coef)
  expect_equal(b$oof$mu[te], pred)
  # fits run under their own seeds: point seed, fold seed + 1000 + k, boot seed + b
  expect_equal(b$point$draw_probe, withr::with_seed(7, stats::runif(1)))
  expect_equal(b$folds[[k]]$draw_probe, withr::with_seed(7 + 1000L + k, stats::runif(1)))
  idx <- withr::with_seed(7 + 1, {
    cc <- train$country
    out <- unlist(lapply(sort(unique(cc)), function(g) { w <- which(cc == g); w[sample.int(length(w), length(w), replace = TRUE)] }))
    out
  })
  expect_equal(b$boot[[1]]$ids, train$ID[idx])
  expect_equal(b$boot[[1]]$draw_probe, withr::with_seed(7 + 1, {
    cc <- train$country
    for (g in sort(unique(cc))) { w <- which(cc == g); sample.int(length(w), length(w), replace = TRUE) }
    stats::runif(1)
  }))
})

test_that("bootstrap keeps the country counts", {
  train <- make_synthetic_data(200, 1, countries = c("DE", "UK"))
  b <- .fit_bundle(toy_method(), train, "S4", "UK", list(), 3, 3, 1)
  for (f in b$boot) {
    expect_equal(table(substr(f$ids, 1, 2))[c("DE", "UK")], table(train$country)[c("DE", "UK")])
  }
})

test_that("the same seed gives the same bundle, and toy bundles save silently", {
  train <- make_synthetic_data(150, 3)
  b1 <- .fit_bundle(toy_method(), train, "S4", "UK", list(), 3, 3, 5)
  b2 <- .fit_bundle(toy_method(), train, "S4", "UK", list(), 3, 3, 5)
  expect_identical(b1$boot, b2$boot); expect_identical(b1$oof, b2$oof)
  expect_identical(b1$folds, b2$folds)
  expect_silent(saveRDS(b1, tempfile()))
})

test_that("uses_oof: point gets the pairs, folds get FALSE, replicates get oof[idx, ]", {
  train <- make_synthetic_data(200, 1)
  b <- .fit_bundle(toy_method(uses_oof = TRUE), train, "S4", "UK", list(), 3, 2, 7)
  o <- b$point$oof_given
  expect_s3_class(o, "data.frame"); expect_equal(nrow(o), 200L)
  expect_named(o, c("oof_mu", "oof_resid"))
  y <- .score(as.matrix(train[DIMS]), b$vs)
  expect_equal(o$oof_mu, b$oof$mu)
  expect_equal(o$oof_resid, y - b$oof$mu)
  for (f in b$folds) expect_identical(f$oof_given, FALSE)
  for (f in b$boot) {
    idx <- match(f$ids, train$ID)
    expect_equal(f$oof_given, o[idx, ], ignore_attr = TRUE)
  }
})

toy_resp_method <- function() list(
  id = "toyr", family = "response", sets = "S4", variants = list(list()),
  aliases = "toyr", needs = character(0), predict_needs = character(0),
  fit = function(train, set, variant, oof = NULL) {
    stopifnot(is.null(train$y))
    list(set = set, variant = variant, n_train = nrow(train))
  },
  simulate = function(fit, newdata, M, key) {
    array(sample(1:5, nrow(newdata) * M * 5, replace = TRUE), c(nrow(newdata), M, 5))
  }
)

test_that("response methods give oof profiles and no target", {
  toy_resp <- toy_resp_method()
  b <- .fit_bundle(toy_resp, make_synthetic_data(120, 2), "S4", NA, list(), 3, 2, 1, Q = 7)
  expect_true(is.na(b$target)); expect_null(b$vs)
  expect_null(b$oof$profiles)
  expect_equal(dim(b$oof$state), c(120L, 7L)); expect_true(is.integer(b$oof$state))
  expect_equal(sum(b$state_freq), 120)
})

test_that("the registry resolves ids and aliases and lists methods on error", {
  m <- toy_method()
  .register_p2e_method(m)
  withr::defer(rm(list = "toy", envir = .p2e_methods))
  expect_equal(.p2e_method("toy")$id, "toy")
  expect_error(.p2e_method("nope"), "toy")
})

test_that("a stored bundle loads from fit_dir and cache_dir", {
  b <- .fit_bundle(toy_method(), make_synthetic_data(60, 1), "S4", "UK", list(), 3, 2, 1)
  d <- tempfile(); dir.create(d); dir.create(file.path(d, "fits"))
  saveRDS(b, file.path(d, "fits", "toy_S4_UK.rds"))
  withr::local_options(promis2eq5d.fit_dir = NULL, promis2eq5d.cache_dir = d)
  expect_s3_class(.load_bundle("toy", "S4", "UK"), "p2e_bundle")
  expect_equal(.cache_dir(), d)
})

test_that("non-converged or failing point and fold fits stop the bundle", {
  tr <- make_synthetic_data(120, 1)
  # call 1 = point fit, call 1 + k = fold k
  expect_error(.fit_bundle(toy_method(bad_call = 3), tr, "S4", "UK", list(), 3, 0, 1),
               "toy.*S4.*fold 2 fit did not converge")
  expect_error(.fit_bundle(toy_method(err_call = 1), tr, "S4", "UK", list(), 3, 0, 1),
               "point fit failed.*errored on purpose")
  expect_error(.fit_bundle(toy_method(bad_call = 1), tr, "S4", "UK", list(), 3, 0, 1),
               "point fit did not converge")
  expect_error(.fit_bundle(toy_method(err_call = 3), tr, "S4", "UK", list(), 3, 0, 1),
               "fold 2 fit failed")
  # uses_oof: the folds come first, the point fit is call K + 1
  expect_error(.fit_bundle(toy_method(uses_oof = TRUE, bad_call = 4), tr, "S4", "UK", list(), 3, 0, 1),
               "point fit did not converge")
  expect_error(.fit_bundle(toy_method(uses_oof = TRUE, bad_call = 2), tr, "S4", "UK", list(), 3, 0, 1),
               "fold 2 fit did not converge")
})
test_that("the shipped methods resolve through the registry", {
  for (id in c("R06_mvop_dim", "R09_chain", "D03_clad", "D12_klapproth")) expect_identical(.p2e_method(id)$id, id)
})

test_that("fit_info$converged = FALSE (no top-level converged) is honoured everywhere", {
  tr <- make_synthetic_data(120, 1)
  expect_error(.fit_bundle(toy_method(fi_bad_call = 1), tr, "S4", "UK", list(), 3, 0, 1),
               "point fit did not converge")
  expect_error(.fit_bundle(toy_method(fi_bad_call = 3), tr, "S4", "UK", list(), 3, 0, 1),
               "fold 2 fit did not converge")
  expect_error(.fit_bundle(toy_method(uses_oof = TRUE, fi_bad_call = 2), tr, "S4", "UK", list(), 3, 0, 1),
               "fold 2 fit did not converge")
  b <- .fit_bundle(toy_method(fi_bad_call = 6), tr, "S4", "UK", list(), 3, 3, 1)  # boot 2 = call 1 + K + 2
  expect_equal(b$boot_failed, 1L); expect_length(b$boot, 3)
  expect_true(.fit_converged(list()))
  expect_false(.fit_converged(list(fit_info = list(converged = NA))))
})

test_that("oof profiles of a response bundle are seeded and leave the caller's stream alone", {
  tr <- make_synthetic_data(120, 2)
  fit <- function() .fit_bundle(toy_resp_method(), tr, "S4", NA, list(), 3, 0, 4, Q = 6)
  set.seed(1); b1 <- fit()
  set.seed(2); before <- .Random.seed; b2 <- fit()
  expect_identical(b1$oof$state, b2$oof$state)
  expect_identical(.Random.seed, before)
  # fold k is simulated under seed + 2000 + k
  te <- which(b1$fold_id == 2L)
  hand <- withr::with_seed(4 + 2000L + 2L, toy_resp_method()$simulate(NULL, tr[te, ], 6, seq_along(te)))
  expect_identical(b1$oof$state[te, ], .state_index(hand))
})

test_that("fit_bundle validates the folds and the bootstrap strata", {
  tr <- make_synthetic_data(120, 1)
  fit <- function(train, ...) .fit_bundle(toy_method(), train, "S4", "UK", list(), 3, 2, 1, ...)
  nostrata <- tr[setdiff(names(tr), c("country", "util_band"))]
  expect_error(make_cv_folds(nostrata, v = 3, seed = 1), "country")
  expect_error(fit(nostrata), "country")
  # folds that are valid in range but leave a fold empty
  expect_error(fit(tr, fold_id = rep(1:2, length.out = 120)), "empty")
  # too few rows for K folds
  expect_error(fit(tr[1:2, ]), "util_band")
  na_country <- tr; na_country$country[3] <- NA
  expect_error(fit(na_country), "country.*missing")
  expect_error(fit(na_country, fold_id = rep(1:3, length.out = 120)), "country.*missing")
})

test_that("fit_bundle requires a country column even with util_band", {
  tr <- make_synthetic_data(120, 1)
  nc <- tr[setdiff(names(tr), "country")]
  expect_error(.fit_bundle(toy_method(), nc, "S4", "UK", list(), 3, 2, 1,
                           fold_id = rep(1:3, length.out = 120)), "country column")
})

# ---- chunked bootstrap: core + attempts in chunks + assemble == .fit_bundle() ----
# meta$created is the only field that depends on when a bundle is assembled
no_created <- function(b) { b$meta$created <- NULL; b }
# a response toy whose replicates depend on their sample and seed
toy_resp_probe <- function() {
  m <- toy_resp_method()
  m$fit <- function(train, set, variant, oof = NULL) {
    stopifnot(is.null(train$y))
    list(set = set, ids = train$ID, folds_seen = train$fold, probe = stats::runif(1))
  }
  m
}
# fails replicate b of a bundle (seed, train) whatever the order or task it is fitted in
fail_sample <- function(m, train, seed, b) {
  bad <- train$ID[.with_fit_seed(seed + b, .boot_index(as.character(train$country), nrow(train)))]
  f <- m$fit
  m$fit <- function(train, ...) { if (identical(train$ID, bad)) stop("toy fit failed on purpose"); f(train, ...) }
  m
}

test_that("core + attempts in scrambled chunks + assemble is identical to .fit_bundle()", {
  withr::local_options(promis2eq5d.repo_commit = "abc", promis2eq5d.repo_dirty = FALSE)
  train <- make_synthetic_data(200, 1)
  fid <- make_cv_folds(train, v = 3, seed = 9)$fold
  B <- 7L
  cfgs <- list(
    direct = list(m = toy_method, target = "UK", fold_id = fid, Q = 50),
    response = list(m = toy_resp_probe, target = NA, fold_id = fid, Q = 4),
    uses_oof = list(m = function() toy_method(uses_oof = TRUE), target = "UK", fold_id = NULL, Q = 50))
  for (nm in names(cfgs)) {
    cf <- cfgs[[nm]]
    ref <- .fit_bundle(cf$m(), train, "S4", cf$target, list(), K = 3, B = B, seed = 7, Q = cf$Q,
                       fold_id = cf$fold_id)
    uses_oof <- nm == "uses_oof"
    chunk <- function(b, core = NULL) {
      .boot_attempts(cf$m(), train, "S4", cf$target, list(), seed = 7, b = b,
                     fold_id = if (is.null(core)) cf$fold_id else core$fold_id,
                     vs = core$vs, pool = core$pool)
    }
    # attempts that need no pool can be made before the core
    c3 <- if (!uses_oof) chunk(c(7L, 6L))
    core <- .bundle_core(cf$m(), train, "S4", cf$target, list(), K = 3, B = B, seed = 7, Q = cf$Q,
                         fold_id = cf$fold_id)
    expect_s3_class(core, "p2e_bundle_core")
    if (uses_oof) c3 <- chunk(c(7L, 6L), core)
    c1 <- chunk(1:3, core); c2 <- chunk(c(5L, 4L), core)
    expect_length(c2, 2L); expect_identical(c2[[1]]$b, 5L)
    got <- .bundle_assemble(cf$m(), train, "S4", cf$target, list(), core, c(c3, c1, c2), B = B, seed = 7)
    expect_identical(no_created(got), no_created(ref), label = nm)
    expect_s3_class(got$meta$created, "POSIXct")
  }
})

test_that("a failed attempt in a chunk is replaced at assembly exactly as in one task", {
  train <- make_synthetic_data(200, 1)
  ref <- .fit_bundle(toy_method(fail_on = 2), train, "S4", "UK", list(), 3, 4, 1)
  # one toy_method(fail_on = 2) through core, chunks and assembly (it counts its calls)
  m <- toy_method(fail_on = 2)
  core <- .bundle_core(m, train, "S4", "UK", list(), 3, 4, 1)
  a12 <- .boot_attempts(m, train, "S4", "UK", list(), 1, 1:2, core$fold_id, core$vs)
  a3 <- .boot_attempts(m, train, "S4", "UK", list(), 1, 3L, core$fold_id, core$vs)
  a4 <- .boot_attempts(m, train, "S4", "UK", list(), 1, 4L, core$fold_id, core$vs)
  expect_null(a12[[2]]$fit); expect_match(a12[[2]]$reason, "failed on purpose")
  got <- .bundle_assemble(m, train, "S4", "UK", list(), core, c(a4, a12, a3), 4, 1)
  expect_identical(no_created(got), no_created(ref))
  expect_equal(got$boot_attempted, 5L)
  # a failure keyed on the sample, every piece with a fresh method, made in scrambled order
  for (uo in c(FALSE, TRUE)) {
    mk <- function() fail_sample(toy_method(uses_oof = uo), train, 1, 3L)
    ref <- .fit_bundle(mk(), train, "S4", "UK", list(), 3, 4, 1)
    core <- .bundle_core(mk(), train, "S4", "UK", list(), 3, 4, 1)
    a <- c(.boot_attempts(mk(), train, "S4", "UK", list(), 1, c(4L, 3L), core$fold_id, core$vs, core$pool),
           .boot_attempts(mk(), train, "S4", "UK", list(), 1, 1:2, core$fold_id, core$vs, core$pool))
    got <- .bundle_assemble(mk(), train, "S4", "UK", list(), core, a, 4, 1)
    expect_identical(no_created(got), no_created(ref))
    expect_identical(got$boot_attempted, ref$boot_attempted)
    expect_identical(got$boot_failed, ref$boot_failed)
    expect_identical(got$boot_failures, ref$boot_failures)
    expect_equal(got$boot_failures$b, 3L)
  }
})

test_that("assembly raises the attempt-cap warning and the all-failed error as in one task", {
  train <- make_synthetic_data(200, 1)
  bad <- function(m, bs) { for (b in bs) m <- fail_sample(m, train, 1, b); m }
  mk <- function() bad(toy_method(), c(2L, 4L, 6L, 7L, 8L))
  core <- .bundle_core(mk(), train, "S4", "UK", list(), 3, 4, 1)
  a <- .boot_attempts(mk(), train, "S4", "UK", list(), 1, 4:1, core$fold_id, core$vs)
  expect_warning(got <- .bundle_assemble(mk(), train, "S4", "UK", list(), core, a, 4, 1), "3 of 4 .*8 attempts")
  expect_warning(ref <- .fit_bundle(mk(), train, "S4", "UK", list(), 3, 4, 1), "3 of 4 .*8 attempts")
  expect_identical(no_created(got), no_created(ref))
  mk <- function() bad(toy_method(), 1:4)
  core <- .bundle_core(mk(), train, "S4", "UK", list(), 3, 2, 1)
  a <- .boot_attempts(mk(), train, "S4", "UK", list(), 1, 1:2, core$fold_id, core$vs)
  expect_error(.bundle_assemble(mk(), train, "S4", "UK", list(), core, a, 2, 1),
               "all 4 bootstrap attempts failed for toy")
})

test_that("assembly needs every attempt 1..B exactly once and the core's settings", {
  train <- make_synthetic_data(200, 1)
  core <- .bundle_core(toy_method(), train, "S4", "UK", list(), 3, 3, 1)
  att <- function(b) .boot_attempts(toy_method(), train, "S4", "UK", list(), 1, b, core$fold_id, core$vs)
  asm <- function(a, B = 3, seed = 1, tr = train) {
    .bundle_assemble(toy_method(), tr, "S4", "UK", list(), core, a, B, seed)
  }
  expect_error(asm(att(c(1L, 3L))), "missing.*\\b2\\b")
  expect_error(asm(c(att(1:3), att(2L))), "duplicated.*\\b2\\b")
  expect_error(asm(att(1:4)), "outside 1\\.\\.3.*\\b4\\b")
  expect_error(asm(att(1:3), B = 4), "B")
  expect_error(asm(att(1:3), seed = 2), "seed")
  expect_error(asm(att(1:3), tr = make_synthetic_data(200, 2)), "training data")
  expect_error(.bundle_assemble(toy_method(), train, "S4", "DE", list(), core, att(1:3), 3, 1), "target UK, not DE")
  expect_error(.bundle_assemble(toy_method(), train, "S4", "UK", list(a = 1), core, att(1:3), 3, 1), "variant")
  rcore <- .bundle_core(toy_resp_probe(), train, "S4", NA, list(), 3, 1, 1, Q = 3)
  ratt <- .boot_attempts(toy_resp_probe(), train, "S4", NA, list(), 1, 1L, rcore$fold_id)
  expect_s3_class(.bundle_assemble(toy_resp_probe(), train, "S4", NA, list(), rcore, ratt, 1, 1), "p2e_bundle")
  expect_error(.bundle_assemble(toy_resp_probe(), train, "S4", "UK", list(), rcore, ratt, 1, 1), "target NA")
  core0 <- .bundle_core(toy_method(), train, "S4", "UK", list(), 3, 0, 1)
  b0 <- .bundle_assemble(toy_method(), train, "S4", "UK", list(), core0, list(), 0, 1)
  expect_length(b0$boot, 0); expect_identical(b0$boot_attempted, 0L)
  expect_error(.bundle_assemble(toy_method(), train, "S4", "UK", list(), core0, att(1L), 0, 1), "outside")
})

test_that(".boot_attempts: a uses_oof method needs the pool; attempts are validated", {
  train <- make_synthetic_data(200, 1)
  fid <- rep(1:3, length.out = 200)
  expect_error(.boot_attempts(toy_method(uses_oof = TRUE), train, "S4", "UK", list(), 1, 1:2, fid),
               "pool")
  expect_error(.boot_attempts(toy_method(), train, "S4", "UK", list(), 1, 0L, fid), "attempt")
  expect_error(.boot_attempts(toy_method(), train, "S4", "UK", list(), 1, 1000L, fid), "attempt")
  expect_error(.boot_attempts(toy_method(), train, "S4", "UK", list(), 1, c(1L, 1L), fid), "attempt")
  expect_error(.boot_attempts(toy_method(), train, "S4", "UK", list(), 1, 1L, fid[-1]), "fold_id")
  expect_identical(.boot_attempts(toy_method(), train, "S4", "UK", list(), 1, integer(), fid), list())
  r <- .boot_attempts(toy_method(), train, "S4", "UK", list(), 1, 2L, fid)
  expect_named(r[[1]], c("b", "fit", "reason")); expect_null(r[[1]]$reason)
  # no attempts need no pool: a uses_oof bundle without folds and B = 0 (an empty pool) still fits
  expect_identical(.boot_attempts(toy_method(uses_oof = TRUE), train, "S4", "UK", list(), 1, integer(), fid), list())
  b0 <- .fit_bundle(toy_method(uses_oof = TRUE), train, "S4", "UK", list(), 3, 0, 1, fit_folds = FALSE)
  expect_s3_class(b0, "p2e_bundle"); expect_equal(nrow(b0$pool), 0L); expect_length(b0$boot, 0)
})

test_that("the core holds no training data frame and no T-score column", {
  train <- make_synthetic_data(200, 1)
  train$T_PF <- 50; train$T_ANX <- 50   # stand-ins for T-score columns (named as .tscores_internal names them)
  for (uo in c(FALSE, TRUE)) {
    core <- .bundle_core(toy_method(uses_oof = uo), train, "S4", "UK", list(), 3, 2, 1)
    expect_false(any(vapply(core, is.data.frame, NA) & names(core) != "pool"))
    # nor one level down (the toy's point fit keeps the pool it was given: that is the pool)
    nested <- unlist(lapply(setdiff(names(core), "pool"), function(k) {
      e <- core[[k]]
      if (!is.list(e) || is.data.frame(e)) return(character())
      hit <- vapply(e, function(x) is.data.frame(x) && !identical(x, core$pool), NA)
      if (any(hit)) paste(k, which(hit)) else character()
    }))
    expect_identical(nested, character())
    expect_false("train" %in% names(core))
    if (uo) expect_named(core$pool, c("oof_mu", "oof_resid"))
    expect_false(any(grepl("T_PF|T_ANX", names(unlist(core)))))
  }
})
