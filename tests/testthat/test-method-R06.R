test_that("R06 keyed simulation with key = rows has the right shape and level range", {
  b <- stand_in_bundle("R06_mvop_dim", "S3"); nd <- make_synthetic_data(5, 2)
  m <- .p2e_method("R06")
  a <- withr::with_seed(3, m$simulate(b$point, nd, 7, key = 1:5))
  expect_equal(dim(a), c(5L, 7L, 5L)); expect_true(all(a %in% 1:5))
})
test_that("R06 common key gives identical profiles for rows of one person with identical inputs", {
  b <- stand_in_bundle("R06_mvop_dim", "S3"); nd <- make_synthetic_data(2, 2)[c(1, 1), ]
  a <- withr::with_seed(3, .p2e_method("R06")$simulate(b$point, nd, 50, key = c(1L, 1L)))
  expect_identical(a[1, , ], a[2, , ])
})
test_that("R06 mean utility is exact-additive plus remainder for any value set", {
  b <- stand_in_bundle("R06_mvop_dim", "S3"); nd <- make_synthetic_data(20, 2)
  for (cc in c("UK", "USA")) {
    vs <- p2e_value_set(cc)
    mc <- rowMeans(.score(withr::with_seed(9, .p2e_method("R06")$simulate(b$point, nd, 20000, 1:20)), vs))
    expect_equal(.p2e_method("R06")$mean_utility(b$point, nd, vs), mc, tolerance = 0.01)
  }
})
test_that("R06 is registered with the interface fields and its fit carries set and variant", {
  m <- .p2e_method("mvop_dim")
  expect_identical(m$id, "R06_mvop_dim"); expect_identical(m$family, "response")
  expect_identical(m$sets, c("S1", "S3", "S2", "S4"))
  expect_identical(m$needs, c("ordinal", "pbivnorm", "Matrix"))
  b <- stand_in_bundle("R06_mvop_dim", "S3")
  expect_identical(b$point$set, "S3"); expect_identical(b$point$variant, list(form = "conceptual"))
  mg <- m$marginals(b$point, make_synthetic_data(4, 2))
  expect_length(mg, 5L); expect_equal(unname(rowSums(mg[[1]])), rep(1, 4))
})
test_that("R06 mean utility is deterministic and leaves the caller's RNG stream alone", {
  b <- stand_in_bundle("R06_mvop_dim", "S3"); nd <- make_synthetic_data(6, 2)
  vs <- p2e_value_set("USA"); m <- .p2e_method("R06")
  set.seed(11); s0 <- .Random.seed
  u1 <- m$mean_utility(b$point, nd, vs)
  expect_identical(.Random.seed, s0)
  expect_identical(u1, m$mean_utility(b$point, nd, vs))
})
test_that("R06 keyed simulation with key = seq_len(n) uses the pipeline's normals layout", {
  b <- stand_in_bundle("R06_mvop_dim", "S3"); nd <- make_synthetic_data(4, 2); M <- 5
  a <- withr::with_seed(5, .p2e_method("R06")$simulate(b$point, nd, M, key = 1:4))
  Z <- withr::with_seed(5, matrix(rnorm(4 * M * 5), 4 * M, 5) %*% chol(b$point$R))
  cm <- .mvop_cuts(b$point, nd)
  lvl1 <- matrix(1L + (Z[, 1] > cm[[1]][rep(1:4, M), 1]), 4, M)
  expect_identical(a[, , 1] >= 2L, lvl1 >= 2L)
})
test_that("a saved R06 bundle does not drag the training data along", {
  b <- stand_in_bundle("R06_mvop_dim", "S3")
  f <- tempfile(fileext = ".rds"); saveRDS(b, f); on.exit(unlink(f))
  message("R06 stand-in bundle size: ", file.size(f), " bytes")
  expect_lt(file.size(f), 5e6)
})
