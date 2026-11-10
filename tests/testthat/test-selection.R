# Variant selection (R/selection.R) against the frozen pipeline
# (fixtures/pipeline_selection.rds, data-raw/pipeline_fixtures.R) and on toy methods.

fx <- readRDS(test_path("fixtures", "pipeline_selection.rds"))

test_that("the R09 grid is the pipeline's 30 configurations", {
  expect_identical(.r09_configs(), fx$ml_configs); expect_length(.r09_configs(), 30L)
})

test_that("the R09 grid leaves the caller's RNG alone and is R09's variant list, default first", {
  set.seed(4); s0 <- .Random.seed
  g <- .r09_configs()
  expect_identical(.Random.seed, s0)
  m <- .p2e_method("R09")
  expect_true(is.function(m$variants))
  expect_identical(.p2e_method_variants(m, "S3"), g)
  expect_identical(.p2e_method_variants(m, "S1")[[1]], g[[1]])
})

test_that("the selection rule matches the pipeline on three tables", {
  for (i in seq_along(fx$tables)) {
    s <- .select_variant(fx$tables[[i]])
    expect_identical(s$variant_id, fx$selected[[i]])
    expect_identical(as.integer(s$n_tied), fx$n_tied[[i]])
    expect_identical(as.integer(s$n_candidates), fx$n_candidates[[i]])
  }
})

test_that("failed variants are never chosen, and a cell with none usable is an error", {
  u <- data.frame(variant_id = 1:3, ok = c(FALSE, TRUE, TRUE), mae = c(0.01, 0.05, 0.06), tie_metric = 1:3)
  expect_identical(.select_variant(u)$variant_id, 2L)
  u$ok <- FALSE; expect_error(.select_variant(u), "no usable variant")
  expect_error(.select_variant(u, "D03_clad S3 UK"), "no usable variant.*D03_clad S3 UK")
})

test_that("a missing MAE among the usable variants is an error; a failed variant's NA is not", {
  u <- data.frame(variant_id = 1:3, ok = TRUE, mae = c(0.05, NA, 0.06), tie_metric = 1:3)
  expect_error(.select_variant(u), "MAE")
  u$ok[2] <- FALSE
  expect_identical(.select_variant(u)$variant_id, 1L)
})

test_that("the tie-break goes tie_metric, then MAE, then variant_id, whatever the row order", {
  u <- data.frame(variant_id = c(3L, 1L, 2L), ok = TRUE, mae = c(0.1004, 0.1002, 0.1000), tie_metric = c(1, 1, 2))
  expect_identical(.select_variant(u)$variant_id, 1L)          # 3 and 1 share tie_metric; 1 has the lower MAE
  expect_identical(.select_variant(u[3:1, ])$variant_id, 1L)
  u$mae[1:2] <- 0.1002
  expect_identical(.select_variant(u)$variant_id, 1L)          # same MAE: the lower variant_id
})

test_that("CRPS and R09 log-probabilities match the pipeline", {
  expect_equal(.crps_sample(fx$crps_y, fx$crps_draws), fx$crps, tolerance = 1e-12)
  skip_if_not_installed("xgboost")
  r09 <- readRDS(test_path("fixtures", "pipeline_R09.rds"))
  m <- .p2e_method("R09")
  fit <- withr::with_seed(r09$fit_seed, m$fit(r09$train, r09$set, r09$variant, NULL))
  nd <- r09$newdata[fx$r09_rows, , drop = FALSE]
  expect_equal(m$logp(fit, nd), log(fx$r09_prob), tolerance = 1e-10)
  expect_equal(m$logp(fit, nd), fx$r09_logp, tolerance = 1e-10)
})

test_that("CRPS handles one draw per row and one row, and refuses missing draws", {
  expect_equal(.crps_sample(c(0.2, 0.5), matrix(c(0.5, 0.5), 2)), c(0.3, 0))
  d <- c(0.1, 0.6, 0.9)
  expect_equal(.crps_sample(0.5, matrix(d, 1)), mean(abs(d - 0.5)) - mean(abs(outer(d, d, "-"))) / 2)
  expect_error(.crps_sample(0.5, matrix(c(0.1, NA), 1)), "missing")
})

test_that("a cv unit reports a fold failure as not usable", {
  tr <- make_synthetic_data(90, 3); fid <- rep(1:3, 30)
  u <- .cv_unit(toy_method(err_call = 2), tr, "S4", "UK", list(), fid, seed = 1)
  expect_false(u$ok)
  expect_match(u$reason, "fold 2")
  expect_match(u$reason, "toy fit errored on purpose")
  expect_true(is.na(u$mae)); expect_true(is.na(u$tie_metric))
  u2 <- .cv_unit(toy_method(bad_call = 3), tr, "S4", "UK", list(), fid, seed = 1)
  expect_false(u2$ok); expect_match(u2$reason, "fold 3.*did not converge")
})

test_that("a direct cv unit uses the bundle's fold fits: MAE of the clamped OOF points, CRPS of seeded draws", {
  tr <- make_synthetic_data(120, 3); fid <- rep(1:3, 40); vs <- p2e_value_set("UK")
  u <- .cv_unit(toy_method(), tr, "S4", "UK", list(), fid, seed = 5)
  b <- .fit_bundle(toy_method(), tr, "S4", "UK", list(), K = 3, B = 0, seed = 5, fold_id = fid)
  y <- .score(as.matrix(tr[DIMS]), vs)
  expect_true(u$ok); expect_true(is.na(u$reason)); expect_identical(u$n, 120L)
  expect_equal(u$mae, mean(abs(y - .clamp(b$oof$mu, vs))), tolerance = 1e-14)
  crps <- numeric(120)
  for (k in 1:3) {
    te <- fid == k
    d <- .with_fit_seed(5 + 2000 + k, toy_method()$draw(b$folds[[k]], tr[te, ], M_DRAWS_SELECT, seq_len(sum(te)),
                                                        "as_published", .vs_min(vs)))
    crps[te] <- .crps_sample(y[te], .clamp(d, vs))
  }
  expect_equal(u$tie_metric, mean(crps), tolerance = 1e-14)
  expect_identical(.cv_unit(toy_method(), tr, "S4", "UK", list(), fid, seed = 5), u)
})

test_that("cv fold fits are the bundle's fold fits (same seed and rows), oof = NULL for uses_oof methods", {
  tr <- make_synthetic_data(120, 3); fid <- rep(1:3, 40)
  tr$y <- .score(as.matrix(tr[DIMS]), p2e_value_set("UK"))
  b <- .fit_bundle(toy_method(uses_oof = TRUE), tr, "S4", "UK", list(), K = 3, B = 0, seed = 5, fold_id = fid)
  for (k in 1:3) {
    r <- .cv_fold_fit(toy_method(uses_oof = TRUE), transform(tr, fold = fid), "S4", list(), fid, k, seed = 5)
    expect_identical(r$fit$coef, b$folds[[k]]$coef)
    expect_identical(r$fit$draw_probe, b$folds[[k]]$draw_probe)
    expect_identical(r$fit$ids, b$folds[[k]]$ids)
    expect_null(r$fit$oof_given)                                # the fold fit computes its own pool
  }
})

test_that("D03 cv fold fits share the bundle fold fit's coefficients and carry their own nested pool", {
  skip_if_not_installed("quantreg")
  tr <- make_synthetic_data(240, 3); m <- .p2e_method("D03")
  tr$util_band <- .util_band(as.matrix(tr[DIMS]), tr$country)
  fid <- make_cv_folds(tr, v = 3, seed = 2)$fold
  b <- .fit_bundle(m, tr, "S3", "UK", list(form = "linear"), K = 3, B = 0, seed = 2)
  expect_identical(b$fold_id, fid)
  tr$y <- .score(as.matrix(tr[DIMS]), p2e_value_set("UK")); tr$fold <- fid
  r <- .cv_fold_fit(m, tr, "S3", list(form = "linear"), fid, 2L, seed = 2)
  expect_equal(r$fit$beta, b$folds[[2]]$beta, tolerance = 1e-14)
  expect_length(r$fit$oof_mu, sum(fid != 2)); expect_equal(sort(unique(r$fit$foldid)), c(1, 3))
  u <- .cv_unit(m, tr, "S3", "UK", list(form = "linear"), fid, seed = 2)
  expect_true(u$ok); expect_true(is.finite(u$mae) && is.finite(u$tie_metric))
})

test_that("a response cv unit scores the OOF mean on each row's own value set and uses -logp", {
  toyr <- list(id = "toyr", family = "response", sets = "S4", variants = list(list()), aliases = "toyr",
               needs = character(0), predict_needs = character(0),
               fit = function(train, set, variant, oof = NULL) list(n = nrow(train)),
               mean_utility = function(fit, newdata, vs) rep(if (vs$code == "UK") 0.8 else 0.7, nrow(newdata)),
               logp = function(fit, newdata) rep(log(0.25), nrow(newdata)))
  tr <- make_synthetic_data(60, 3, countries = c("DE", "UK")); fid <- rep(1:3, 20)
  u <- .cv_unit(toyr, tr, "S4", NA, list(), fid, seed = 1)
  own <- vapply(seq_len(60), function(i) .score(as.matrix(tr[i, DIMS]), p2e_value_set(tr$country[i])), 0)
  expect_true(u$ok)
  expect_equal(u$mae, mean(abs(own - ifelse(tr$country == "UK", 0.8, 0.7))), tolerance = 1e-14)
  expect_equal(u$tie_metric, -log(0.25))
  toyr$logp <- NULL
  expect_true(is.na(.cv_unit(toyr, tr, "S4", NA, list(), fid, seed = 1)$tie_metric))
})

test_that("cv unit checks fold_id", {
  tr <- make_synthetic_data(30, 3)
  expect_error(.cv_unit(toy_method(), tr, "S4", "UK", list(), rep(1:3, 9), seed = 1), "fold_id")
  expect_error(.cv_unit(toy_method(), tr, "S4", "UK", list(), rep(c(1L, 3L), 15), seed = 1), "fold_id")
})
