test_that("D03 point is the capped median and draws come from residual pools", {
  skip_if_not_installed("quantreg")
  bd <- stand_in_bundle("D03_clad", "S3", "UK"); nd <- make_synthetic_data(6, 2)
  m <- .p2e_method("D03")
  expect_true(all(m$point(bd$point, nd) <= 1))
  d <- withr::with_seed(4, m$draw(bd$point, nd, 9, 1:6, "corrected", .vs_min(p2e_value_set("UK")),
                                           pool = bd$pool))
  expect_equal(dim(d), c(6L, 9L))
  same <- withr::with_seed(4, residual_draws(m$point(bd$point, nd), bd$pool$oof_mu, bd$pool$oof_resid, 9))
  expect_identical(d, same)
})
test_that("D03 keyed draws reproduce the documented scheme", {
  skip_if_not_installed("quantreg")
  bd <- stand_in_bundle("D03_clad", "S3", "UK"); m <- .p2e_method("D03")
  nd <- make_synthetic_data(5, 2)[c(1, 1, 2, 3, 3), ]; key <- c(1L, 1L, 2L, 3L, 3L); M <- 20
  d <- withr::with_seed(5, m$draw(bd$point, nd, M, key, "corrected", -0.5, pool = bd$pool))
  expect_equal(dim(d), c(5L, M))
  expect_identical(d[1, ], d[2, ]); expect_identical(d[4, ], d[5, ])
  expect_false(isTRUE(all.equal(d[1, ], d[3, ])))
  mu <- m$point(bd$point, nd)
  U <- withr::with_seed(5, matrix(runif(max(key) * M), max(key), M))
  g <- .rd_groups(mu, bd$pool$oof_mu)
  hand <- matrix(NA_real_, 5, M)
  for (i in 1:5) {
    pool <- sort(bd$pool$oof_resid[g$bin %in% as.integer(strsplit(g$key[i], ",")[[1]])])
    hand[i, ] <- mu[i] + pool[pmax(1, ceiling(U[key[i], ] * length(pool)))]
  }
  expect_equal(d, hand)
})
test_that("D03 interface fields and variants by set", {
  m <- .p2e_method("clad")
  expect_identical(m$id, "D03_clad"); expect_identical(m$family, "direct")
  expect_identical(m$sets, c("S3", "S4", "S1", "S2")); expect_true(m$uses_oof)
  expect_identical(m$needs, c("quantreg", "glmnet")); expect_identical(m$predict_needs, character(0))
  expect_identical(.p2e_method_variants(m, "S1")[[1]], list(form = "numeric"))
  expect_identical(.p2e_method_variants(m, "S2")[[1]], list(form = "numeric"))
  expect_identical(.p2e_method_variants(m, "S3")[[1]], list(form = "linear"))
  expect_identical(.p2e_method_variants(m, "S4")[[1]], list(form = "linear"))
  expect_identical(vapply(.p2e_method_variants(m, "S3"), `[[`, "", "form"),
                   c("linear", "quadratic", "spline3", "spline3_int"))
})
test_that(".p2e_method_variants accepts a list or a function and validates", {
  expect_identical(.p2e_method_variants(list(variants = list(list(a = 1))), "S1"), list(list(a = 1)))
  expect_identical(.p2e_method_variants(list(variants = function(set) list(list(s = set))), "S4"), list(list(s = "S4")))
  expect_error(.p2e_method_variants(list(variants = list()), "S1"), "variants")
  expect_error(.p2e_method_variants(list(variants = function(set) NULL), "S1"), "variants")
})
test_that("D03 stored fit is small and holds no crq object", {
  skip_if_not_installed("quantreg")
  bd <- stand_in_bundle("D03_clad", "S3", "UK")
  f <- file.path(getOption("promis2eq5d.fit_dir"), "D03_clad_S3_UK.rds")
  expect_lt(file.size(f), 5e6)
  expect_identical(bd$point$set, "S3"); expect_identical(bd$point$variant, list(form = "linear"))
  expect_null(bd$point[["fit"]])
})
test_that("bundle OOF equals the standalone fit's internal OOF under the same folds", {
  skip_if_not_installed("quantreg")
  m <- .p2e_method("D03")
  train <- make_synthetic_data(120, seed = 2, countries = c("DE", "UK"))
  b <- .fit_bundle(m, train, "S3", "UK", list(form = "linear"), K = 3, B = 0, seed = 7)
  train$y <- .score(as.matrix(train[DIMS]), p2e_value_set("UK")); train$fold <- b$fold_id
  s <- m$fit(train, "S3", list(form = "linear"), oof = NULL)
  expect_equal(s$oof_mu, b$pool$oof_mu); expect_equal(s$oof_resid, b$pool$oof_resid)
  expect_equal(s$foldid, b$fold_id)
  expect_equal(s$beta, b$point$beta)
})
test_that("D03 with oof = FALSE fits the core model only and cannot draw", {
  skip_if_not_installed("quantreg")
  m <- .p2e_method("D03"); bd <- stand_in_bundle("D03_clad", "S3", "UK")
  fo <- bd$folds[[1]]
  expect_null(fo$oof_mu)
  nd <- make_synthetic_data(4, 2)
  expect_length(m$point(fo, nd), 4L)
  expect_error(m$draw(fo, nd, 5, 1:4, "corrected", -0.5), "cannot draw")
})
test_that("D03 lasso_dummies on S1 fits and predicts", {
  skip_if_not_installed("quantreg"); skip_if_not_installed("glmnet")
  m <- .p2e_method("D03")
  train <- make_synthetic_data(150, seed = 3, countries = c("DE", "UK"))
  train$y <- .score(as.matrix(train[DIMS]), p2e_value_set("UK"))
  train$fold <- make_cv_folds(train, v = 3, seed = 1)$fold
  f <- m$fit(train, "S1", list(form = "lasso_dummies"), oof = FALSE)
  expect_true(isTRUE(f$converged))
  p <- m$point(f, make_synthetic_data(5, 2))
  expect_length(p, 5L); expect_true(all(is.finite(p)) && all(p <= 1))
})
test_that("lasso_dummies bundle OOF equals the standalone fit's internal OOF", {
  skip_if_not_installed("quantreg"); skip_if_not_installed("glmnet")
  m <- .p2e_method("D03"); v <- list(form = "lasso_dummies")
  train <- make_synthetic_data(120, seed = 4, countries = c("DE", "UK"))
  b <- .fit_bundle(m, train, "S1", "UK", v, K = 4, B = 0, seed = 7)
  train$y <- .score(as.matrix(train[DIMS]), p2e_value_set("UK")); train$fold <- b$fold_id
  s <- m$fit(train, "S1", v, oof = NULL)
  expect_equal(s$oof_mu, b$pool$oof_mu); expect_equal(s$oof_resid, b$pool$oof_resid)
  expect_equal(s$cols, b$point$cols); expect_equal(s$beta, b$point$beta)
})

test_that("D03 common key gives positively rank-correlated draws across different pools", {
  skip_if_not_installed("quantreg")
  bd <- stand_in_bundle("D03_clad", "S3", "UK"); m <- .p2e_method("D03")
  nd <- make_synthetic_data(300, 2)
  mu <- m$point(bd$point, nd)
  g <- .rd_groups(mu, bd$pool$oof_mu)
  expect_gt(length(unique(g$key)), 1L)
  # one row from each of two different pools, sharing one key (one person, two time points)
  i <- c(which(g$key == g$key[1])[1], which(g$key != g$key[1])[1])
  M <- 2000
  d <- withr::with_seed(3, m$draw(bd$point, nd[i, ], M, c(1L, 1L), "corrected", -0.5, pool = bd$pool))
  r <- d - mu[i]
  # ties in the pools keep it below 1; the unsorted pools gave about 0
  expect_gt(stats::cor(r[1, ], r[2, ], method = "spearman"), 0.5)
  # independent rows are not correlated
  d2 <- withr::with_seed(3, m$draw(bd$point, nd[i, ], M, 1:2, "corrected", -0.5, pool = bd$pool))
  expect_lt(abs(stats::cor((d2 - mu[i])[1, ], (d2 - mu[i])[2, ], method = "spearman")), 0.2)
})
test_that("rebuilding each replicate's pool equals drawing with the pool handed back to the fit", {
  skip_if_not_installed("quantreg")
  train <- make_synthetic_data(300, 1)
  b <- .fit_bundle(.p2e_method("D03"), train, "S3", "UK", list(form = "linear"), K = 3, B = 2, seed = 5)
  nd <- make_synthetic_data(8, 2); vs <- p2e_value_set("UK"); key <- 1:8
  got <- withr::with_seed(3, .draw_matrix(b, nd, vs, 6, key, "corrected"))
  fits <- lapply(b$boot, function(f) {
    p <- .boot_pool(b, f); f$oof_mu <- p$oof_mu; f$oof_resid <- p$oof_resid
    f$boot_b <- NULL
    f
  })
  b2 <- b; b2$pool <- NULL
  ref <- withr::with_seed(3, .draw_matrix(b2, nd, vs, 6, key, "corrected", fits = fits))
  expect_identical(got, ref)
})
