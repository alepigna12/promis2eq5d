test_that("a bundle has K fold fits, B boot fits and aligned oof", {
  b <- .fit_bundle(toy_method(), train = make_synthetic_data(200, 1), set = "S4", target = "UK",
                   variant = list(), K = 3, B = 4, seed = 1)
  expect_s3_class(b, "p2e_bundle"); expect_length(b$folds, 3); expect_length(b$boot, 4)
  expect_equal(length(b$oof$mu), 200); expect_equal(sort(unique(b$fold_id)), 1:3)
  expect_equal(sum(b$state_freq), 200); expect_equal(dim(b$y_profiles), c(200L, 5L))
  expect_s3_class(b$vs, "p2e_vs"); expect_equal(b$target, "UK")
  expect_named(b$meta, c("data_doi","seed","repo_commit","ihea_tag","r_version","pkg_versions","created"))
})
test_that("failed bootstrap replicates are dropped and counted", {
  b <- .fit_bundle(toy_method(fail_on = 2), make_synthetic_data(200, 1), "S4", "UK", list(), 3, 4, 1)
  expect_length(b$boot, 3); expect_equal(b$boot_failed, 1L); expect_equal(b$boot_attempted, 4L)
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

test_that("response methods give oof profiles and no target", {
  toy_resp <- list(
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
  b <- .fit_bundle(toy_resp, make_synthetic_data(120, 2), "S4", NA, list(), 3, 2, 1, Q = 7)
  expect_true(is.na(b$target)); expect_null(b$vs)
  expect_equal(dim(b$oof$profiles), c(120L, 7L, 5L))
  expect_true(is.integer(b$oof$profiles))
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
  for (id in c("R06_mvop_dim", "R09_chain", "D03_clad")) expect_identical(.p2e_method(id)$id, id)
})
