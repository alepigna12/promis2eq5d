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
test_that("R09 simulate with key = seq_len(n) matches a step-by-step reference using the pipeline's RNG layout", {
  skip_if_not_installed("xgboost")
  b <- stand_in_bundle("R09_chain", "S3"); f <- b$point; nd <- make_synthetic_data(4, 2); n <- 4L; M <- 6L
  a <- withr::with_seed(5, .p2e_method("R09")$simulate(f, nd, M, key = 1:n))
  ref <- withr::with_seed(5, {
    X <- .r09_x(nd, f$set); Xr <- X[rep(seq_len(n), M), , drop = FALSE]
    lv <- matrix(0L, n * M, 5)                     # chain order; row (i, m) = i + n * (m - 1)
    for (k in 1:5) {
      u <- matrix(runif(n * M), n, M)
      p <- .r09_cond(f, k, Xr, lv[, seq_len(k - 1L), drop = FALSE])
      cp <- t(apply(p, 1, cumsum))
      lv[, k] <- 1L + rowSums(as.vector(u) > cp[, 1:4])
    }
    lv
  })
  expect_identical(a, array(as.integer(ref[, match(DIMS, f$order)]), c(n, M, 5)))
})
test_that("R09 handles a single-row newdata", {
  skip_if_not_installed("xgboost")
  b <- stand_in_bundle("R09_chain", "S3"); nd <- make_synthetic_data(2, 2)[1, , drop = FALSE]
  vs <- p2e_value_set("USA"); m <- .p2e_method("R09")
  expect_equal(dim(withr::with_seed(1, m$simulate(b$point, nd, 4, key = 1L))), c(1L, 4L, 5L))
  mg <- m$marginals(b$point, nd)
  expect_length(mg, 5L)
  for (d in mg) { expect_equal(dim(d), c(1L, 5L)); expect_equal(unname(rowSums(d)), 1) }
  u <- m$mean_utility(b$point, nd, vs)
  expect_length(u, 1L); expect_true(is.finite(u) && u >= .vs_min(vs) && u <= 1)
})
test_that("R09 blocked Monte Carlo mean: deterministic, shaped, and equal to the unblocked run within one block", {
  skip_if_not_installed("xgboost")
  b <- stand_in_bundle("R09_chain", "S3"); f <- b$point; nd <- make_synthetic_data(12, 2)
  vs <- p2e_value_set("USA")
  u3 <- .r09_mean_utility(f, nd, vs, block = 5L)
  expect_length(u3, 12L); expect_identical(u3, .r09_mean_utility(f, nd, vs, block = 5L))
  mg3 <- .r09_marginals(f, nd, block = 5L)
  for (d in mg3) { expect_equal(dim(d), c(12L, 5L)); expect_equal(unname(rowSums(d)), rep(1, 12)) }
  expect_identical(.r09_mean_utility(f, nd, vs, block = 12L), .r09_mean_utility(f, nd, vs, block = 500L))
  ref <- withr::with_seed(MC_SEED, .rng_kind = "Mersenne-Twister", .rng_normal_kind = "Inversion",
                          .rng_sample_kind = "Rejection",
                          rowMeans(.score(.r09_simulate(f, nd, N_MC_MEAN, 1:12), vs)))
  expect_identical(.r09_mean_utility(f, nd, vs), ref)
})
test_that("a saved R09 bundle is small, and reloading it simulates identically (raw boosters)", {
  skip_if_not_installed("xgboost")
  b <- stand_in_bundle("R09_chain", "S3"); nd <- make_synthetic_data(5, 2)
  f <- file.path(getOption("promis2eq5d.fit_dir"), .bundle_file("R09_chain", "S3", NA))
  expect_true(file.exists(f))
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
