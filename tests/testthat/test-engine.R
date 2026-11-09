test_that("draw d uses bootstrap fit ((d - 1) mod B) + 1 for every row", {
  b <- stand_in_bundle("R06_mvop_dim", "S3"); nd <- make_synthetic_data(4, 2)
  r <- withr::with_seed(1, .draw_matrix(b, nd, p2e_value_set("UK"), D = 10, key = 1:4, draw_type = "corrected"))
  B <- length(b$boot); expect_equal(r$fit_index, rep_len(seq_len(B), 10))
})
test_that("common key shares individual noise; independent does not", {
  b <- stand_in_bundle("D12_klapproth", "S3K", "UK"); nd <- make_synthetic_data(2, 2)[c(1, 1), ]
  vs <- p2e_value_set("UK")
  com <- withr::with_seed(1, .draw_matrix(b, nd, vs, 200, key = c(1L, 1L), "corrected"))$u
  ind <- withr::with_seed(1, .draw_matrix(b, nd, vs, 200, key = 1:2, "corrected"))$u
  expect_identical(com[1, ], com[2, ]); expect_false(identical(ind[1, ], ind[2, ]))
})
test_that("crosswalked D03 draws land on the target scale with an empty gap", {
  b <- stand_in_bundle("D03_clad", "S3", "UK"); nd <- make_synthetic_data(50, 2); us <- p2e_value_set("USA")
  k <- .fit_crosswalk(p2e_value_set("UK"), us, b$state_freq)
  u <- withr::with_seed(2, .draw_matrix(b, nd, us, 200, 1:50, "corrected", crosswalk = k))$u
  g <- 1 - .vs_gap(us)
  expect_false(any(u > g & u < 1)); expect_true(all(u >= .vs_min(us) & u <= 1))
})
test_that("as_published warns once per call", {
  b <- stand_in_bundle("D12_klapproth", "S3K", "UK"); nd <- make_synthetic_data(3, 2)
  n <- 0
  withCallingHandlers(.draw_matrix(b, nd, p2e_value_set("UK"), 5, 1:3, "as_published"), warning = function(w) {
    n <<- n + 1; expect_match(conditionMessage(w), "as published.*biased"); invokeRestart("muffleWarning") })
  expect_equal(n, 1)
})
test_that("as_published warning names the reason for each method and points to the docs", {
  nd <- make_synthetic_data(3, 2); vs <- p2e_value_set("UK")
  w12 <- capture_warning(.draw_matrix(stand_in_bundle("D12_klapproth", "S3K", "UK"), nd, vs, 3, 1:3, "as_published"))
  expect_match(conditionMessage(w12), "clamp")
  expect_match(conditionMessage(w12), "predict_eq5d")
  w03 <- capture_warning(.draw_matrix(stand_in_bundle("D03_clad", "S3", "UK"), nd, vs, 3, 1:3, "as_published"))
  expect_match(conditionMessage(w03), "gap")
  expect_match(conditionMessage(w03), "predict_eq5d")
})
test_that("one row and a block boundary inside a key both work", {
  b <- stand_in_bundle("R06_mvop_dim", "S3"); nd <- make_synthetic_data(2, 2)[1, , drop = FALSE]
  expect_equal(dim(.draw_matrix(b, nd, p2e_value_set("UK"), 3, 1L, "corrected")$u), c(1L, 3L))
  bd <- stand_in_bundle("D03_clad", "S3", "UK")
  expect_equal(dim(.draw_matrix(bd, nd, p2e_value_set("UK"), 3, 1L, "corrected")$u), c(1L, 3L))
  nd2 <- make_synthetic_data(200, 3); key <- rep(1:200, each = 3)  # 600 rows, block size 500
  r <- withr::with_seed(1, .draw_matrix(b, nd2[key, ], p2e_value_set("UK"), 5, key, "corrected"))
  expect_identical(r$profiles[499, , ], r$profiles[500, , ])       # key 167 spans rows 499-501
  expect_identical(r$profiles[500, , ], r$profiles[501, , ])
})
test_that("fits = list(point) draws every column from the point fit", {
  b <- stand_in_bundle("R06_mvop_dim", "S3"); nd <- make_synthetic_data(4, 2)
  r <- .draw_matrix(b, nd, p2e_value_set("UK"), 6, 1:4, "corrected", fits = list(b$point))
  expect_equal(r$fit_index, rep(1L, 6))
})

test_that("blocks keep non-contiguous keys together and restore row order", {
  testthat::local_mocked_bindings(ENGINE_BLOCK = 4L)
  b <- stand_in_bundle("D12_klapproth", "S3K", "UK"); vs <- p2e_value_set("UK")
  nd0 <- make_synthetic_data(6, 2)
  key <- c(1L, 2L, 3L, 1L, 4L, 2L, 5L, 6L, 3L, 1L)   # 10 rows, keys interleaved
  nd <- nd0[key, ]
  r <- withr::with_seed(3, .draw_matrix(b, nd, vs, 6, key, "corrected"))
  expect_equal(dim(r$u), c(10L, 6L))
  for (k in unique(key)) {
    rows <- which(key == k)
    for (i in rows[-1]) expect_identical(r$u[i, ], r$u[rows[1], ])
  }
  # original order: rows with the same covariates and key agree; distinct keys differ
  expect_false(identical(r$u[1, ], r$u[2, ]))
  expect_equal(.block_rows(key, 4L), list(c(1L, 4L, 10L), c(2L, 3L, 6L, 9L), c(5L, 7L, 8L)))
})
test_that("blocks: a large key forms its own block; keys are numbered by first appearance", {
  expect_equal(.block_rows(c(5L, 5L, 5L, 5L, 5L, 9L), 2L), list(1:5, 6L))
  expect_equal(.block_rows(c(7L, 3L), 500L), list(1:2))
})
test_that("one block equals the unblocked reference (direct, corrected D03, and response)", {
  vs <- p2e_value_set("UK"); nd <- make_synthetic_data(5, 2); key <- 1:5; D <- 7
  bd <- stand_in_bundle("D03_clad", "S3", "UK"); B <- length(bd$boot)
  fi <- rep_len(seq_len(B), D)
  ref <- withr::with_seed(11, {
    u <- matrix(NA_real_, 5, D)
    for (bb in seq_len(B)) {
      cols <- which(fi == bb)
      u[, cols] <- .p2e_method("D03_clad")$draw(bd$boot[[bb]], nd, length(cols), key, "corrected", .vs_min(bd$vs))
    }
    V <- matrix(stats::runif(5 * D), 5, D)[key, , drop = FALSE]
    .clamp(.round_gap(.clamp(u, bd$vs), vs, V), vs)
  })
  got <- withr::with_seed(11, .draw_matrix(bd, nd, vs, D, key, "corrected"))$u
  expect_identical(got, ref)
  br <- stand_in_bundle("R06_mvop_dim", "S3"); Br <- length(br$boot); fir <- rep_len(seq_len(Br), D)
  refp <- withr::with_seed(12, {
    p <- array(NA_integer_, c(5, D, 5))
    for (bb in seq_len(Br)) {
      cols <- which(fir == bb)
      p[, cols, ] <- .p2e_method("R06_mvop_dim")$simulate(br$boot[[bb]], nd, length(cols), key)
    }
    p
  })
  gotr <- withr::with_seed(12, .draw_matrix(br, nd, vs, D, key, "corrected"))
  expect_identical(gotr$profiles, refp)
  expect_identical(gotr$u, .score(refp, vs))
})
test_that("D < B and D not a multiple of B work", {
  b <- stand_in_bundle("D12_klapproth", "S3K", "UK"); vs <- p2e_value_set("UK"); nd <- make_synthetic_data(3, 2)
  B <- length(b$boot)
  r1 <- .draw_matrix(b, nd, vs, 2, 1:3, "corrected")
  expect_equal(dim(r1$u), c(3L, 2L)); expect_equal(r1$fit_index, 1:2); expect_false(anyNA(r1$u))
  r2 <- .draw_matrix(b, nd, vs, B + 1, 1:3, "corrected")
  expect_equal(dim(r2$u), c(3L, B + 1L)); expect_equal(r2$fit_index, c(seq_len(B), 1L)); expect_false(anyNA(r2$u))
  expect_null(r2$profiles)
  br <- stand_in_bundle("R06_mvop_dim", "S3")
  r3 <- .draw_matrix(br, nd, vs, 2, 1:3, "corrected")
  expect_equal(dim(r3$profiles), c(3L, 2L, 5L))
})
test_that("a direct method on another value set without a crosswalk errors", {
  b <- stand_in_bundle("D03_clad", "S3", "UK"); nd <- make_synthetic_data(3, 2)
  expect_error(.draw_matrix(b, nd, p2e_value_set("USA"), 3, 1:3, "corrected"), "crosswalk")
})
test_that("a bad draw_type errors", {
  b <- stand_in_bundle("D03_clad", "S3", "UK"); nd <- make_synthetic_data(3, 2)
  expect_error(.draw_matrix(b, nd, p2e_value_set("UK"), 3, 1:3, "nope"), "draw_type")
  expect_error(.draw_matrix(b, nd, p2e_value_set("UK"), 3, 1:3, NA_character_), "draw_type")
})
