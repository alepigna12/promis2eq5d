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

test_that("oof and bootstrap use the documented folds and indices", {
  train <- make_synthetic_data(200, 1)
  b <- .fit_bundle(toy_method(), train, "S4", "UK", list(), 3, 2, 7)
  expect_equal(b$fold_id, make_cv_folds(train, v = 3, seed = 7)$fold)
  # oof$mu is the prediction of the fold fit that left the row out
  y <- .score(as.matrix(train[DIMS]), b$vs)
  k <- 2L; te <- b$fold_id == k
  expect_equal(b$oof$mu[te], as.numeric(predict(b$folds[[k]]$model, train[te, PREDICTOR_SETS$S4])))
  # bootstrap b resamples rows with sample.int under seed + b
  idx <- withr::with_seed(7 + 1, sample.int(200, 200, replace = TRUE))
  expect_equal(b$boot[[1]]$n_train, 200L)
  tr <- train; tr$y <- y
  expect_equal(unname(coef(b$boot[[1]]$model)),
               unname(coef(lm(y ~ ., data = cbind(tr[idx, PREDICTOR_SETS$S4], y = y[idx])))))
  # no oof for methods that do not ask for it
  expect_null(b$point$oof_given)
})

test_that("uses_oof passes oof to the point fit and resampled oof to replicates", {
  train <- make_synthetic_data(200, 1)
  b <- .fit_bundle(toy_method(uses_oof = TRUE), train, "S4", "UK", list(), 3, 2, 7)
  o <- b$point$oof_given
  expect_named(o, c("oof_mu", "oof_resid"))
  y <- .score(as.matrix(train[DIMS]), b$vs)
  expect_equal(o$oof_mu, b$oof$mu)
  expect_equal(o$oof_resid, y - b$oof$mu)
  expect_null(b$folds[[1]]$oof_given)
  idx <- withr::with_seed(7 + 2, sample.int(200, 200, replace = TRUE))
  expect_equal(b$boot[[2]]$oof_given$oof_mu, o$oof_mu[idx])
  expect_equal(b$boot[[2]]$oof_given$oof_resid, o$oof_resid[idx])
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
  # the toy lm keeps its formula environment, which refers to the loaded package
  suppressWarnings(saveRDS(b, file.path(d, "fits", "toy_S4_UK.rds")))
  withr::local_options(promis2eq5d.fit_dir = NULL, promis2eq5d.cache_dir = d)
  expect_s3_class(.load_bundle("toy", "S4", "UK"), "p2e_bundle")
  expect_equal(.cache_dir(), d)
})
