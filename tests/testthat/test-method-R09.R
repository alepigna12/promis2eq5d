test_that("R09 keyed simulation has the right shape and level range", {
  skip_if_not_installed("xgboost")
  b <- stand_in_bundle("R09_chain", "S3"); nd <- make_synthetic_data(5, 2)
  a <- withr::with_seed(3, .p2e_method("R09")$simulate(b$point, nd, 7, key = 1:5))
  expect_equal(dim(a), c(5L, 7L, 5L)); expect_true(all(a %in% 1:5)); expect_type(a, "integer")
})
test_that("R09 common key gives identical profiles for rows of one person with identical inputs", {
  skip_if_not_installed("xgboost")
  b <- stand_in_bundle("R09_chain", "S3"); nd <- make_synthetic_data(2, 2)[c(1, 1), ]
  a <- withr::with_seed(3, .p2e_method("R09")$simulate(b$point, nd, 50, key = c(1L, 1L)))
  expect_identical(a[1, , ], a[2, , ])
})
test_that("R09 is registered with the interface fields and its fit carries set and variant", {
  skip_if_not_installed("xgboost")
  m <- .p2e_method("chain")
  expect_identical(m$id, "R09_chain"); expect_identical(m$family, "response")
  expect_identical(m$sets, c("S3", "S1", "S4", "S2"))
  expect_identical(m$needs, "xgboost"); expect_identical(m$predict_needs, "xgboost")
  expect_identical(m$variants[[1]], list(eta = 0.3, max_depth = 2, subsample = 1,
                                         colsample_bytree = 1, min_child_weight = 1, nrounds = 10))
  b <- stand_in_bundle("R09_chain", "S3")
  expect_identical(b$point$set, "S3"); expect_identical(b$point$variant, m$variants[[1]])
  expect_identical(b$point$order, c("E4", "E1", "E3", "E2", "E5"))
})
test_that("R09 marginals rows sum to 1", {
  skip_if_not_installed("xgboost")
  b <- stand_in_bundle("R09_chain", "S3"); nd <- make_synthetic_data(4, 2)
  mg <- .p2e_method("R09")$marginals(b$point, nd)
  expect_length(mg, 5L); expect_identical(names(mg), DIMS)
  for (d in mg) { expect_equal(dim(d), c(4L, 5L)); expect_equal(unname(rowSums(d)), rep(1, 4)) }
})
test_that("R09 mean utility lies in [min, 1], is deterministic and leaves the caller's RNG alone", {
  skip_if_not_installed("xgboost")
  b <- stand_in_bundle("R09_chain", "S3"); nd <- make_synthetic_data(6, 2)
  vs <- p2e_value_set("USA"); m <- .p2e_method("R09")
  set.seed(11); s0 <- .Random.seed
  u1 <- m$mean_utility(b$point, nd, vs)
  expect_identical(.Random.seed, s0)
  expect_identical(u1, m$mean_utility(b$point, nd, vs))
  expect_true(all(u1 >= .vs_min(vs) & u1 <= 1))
})
test_that("R09 keyed simulation with key = seq_len(n) draws one n x M block of uniforms per chain step", {
  skip_if_not_installed("xgboost")
  b <- stand_in_bundle("R09_chain", "S3"); nd <- make_synthetic_data(4, 2); M <- 5
  a <- withr::with_seed(5, .p2e_method("R09")$simulate(b$point, nd, M, key = 1:4))
  U1 <- withr::with_seed(5, matrix(runif(4 * M), 4, M))     # first chain step is E4
  X <- .r09_x(nd, "S3")
  p <- .r09_cond(b$point, 1L, X, matrix(0, 4, 0))
  lvl <- matrix(1L, 4, M)
  for (j in 1:4) lvl <- lvl + (U1 > matrix(t(apply(p, 1, cumsum))[, j], 4, M))
  expect_identical(a[, , 4], lvl)
})
test_that("a saved R09 bundle is small, and reloading it simulates identically (raw boosters)", {
  skip_if_not_installed("xgboost")
  b <- stand_in_bundle("R09_chain", "S3"); nd <- make_synthetic_data(5, 2)
  f <- file.path(getOption("promis2eq5d.fit_dir"), "R09_chain_S3_any.rds")
  expect_true(file.exists(f))
  message("R09 stand-in bundle size: ", file.size(f), " bytes")
  expect_lt(file.size(f), 5e6)
  b2 <- readRDS(f); m <- .p2e_method("R09")
  expect_true(is.raw(b2$point$dims$E1$raw))
  a1 <- withr::with_seed(7, m$simulate(b$point, nd, 30, key = 1:5))
  a2 <- withr::with_seed(7, m$simulate(b2$point, nd, 30, key = 1:5))
  expect_identical(a2, a1)
})
test_that("a reloaded booster predicts with one thread", {
  skip_if_not_installed("xgboost")
  b <- stand_in_bundle("R09_chain", "S3")
  bst <- .r09_boosters(b$point)[[1]]
  expect_identical(as.integer(xgboost::xgb.config(bst)$learner$generic_param$nthread), 1L)
})
