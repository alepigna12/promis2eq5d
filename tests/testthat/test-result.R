use_stand_in_fits()
quiet <- function(expr) suppressMessages(expr)

res_data <- function() {
  nd <- user_frame(make_synthetic_data(12, 2)); nd$arm <- rep(c("a", "b"), 6); nd
}

test_that("print shows a header and the first rows (mean mode, R06)", {
  r <- quiet(predict_eq5d(res_data(), "UK"))
  out <- capture.output(print(r))
  expect_match(out[1], "p2e_result")
  expect_true(any(grepl("method: R06_mvop_dim", out)))
  expect_true(any(grepl("set: S1", out)))
  expect_true(any(grepl("value set: UK", out)))
  expect_true(any(grepl("route:", out)))
  expect_true(any(grepl("mode: mean", out)))
  expect_true(any(grepl("utility", out)))
  capture.output(shown <- print(r))
  expect_identical(shown, r)
})

test_that("print shows draws and the crosswalk route (individual, D03)", {
  r <- quiet(predict_eq5d(res_data(), "USA", family = "direct", mode = "individual", draws = 4, seed = 1))
  out <- capture.output(print(r))
  expect_true(any(grepl("draws: 4", out)))
  expect_true(any(grepl("crosswalk from", out)))
  expect_true(any(grepl("method: D03_clad", out)))
  expect_lte(length(out), 20)
})

test_that("summary: mean mode", {
  r <- quiet(predict_eq5d(res_data(), "UK"))
  s <- summary(r)
  expect_s3_class(s, "data.frame")
  expect_named(s, c("n", "mean", "sd", "min", "max"))
  expect_equal(s$n, 12); expect_equal(s$mean, mean(r$utility)); expect_equal(s$max, max(r$utility))
})

test_that("summary: individual mode", {
  nd <- res_data(); nd$pid <- paste0("p", 1:12)
  r <- quiet(predict_eq5d(nd, "UK", mode = "individual", draws = 20, seed = 1, id = "pid"))
  s <- summary(r, level = 0.8)
  expect_named(s, c(".row", "pid", "mean", "sd", "lower", "upper"))
  expect_equal(nrow(s), 12)
  u <- as.matrix(r)
  expect_equal(s$mean, unname(rowMeans(u)))
  expect_equal(s$lower, unname(apply(u, 1, stats::quantile, 0.1)))
  expect_equal(s$upper, unname(apply(u, 1, stats::quantile, 0.9)))
})

test_that("summary: cohort mode uses the 95% level", {
  r <- quiet(predict_eq5d(res_data(), "UK", mode = "cohort", group = "arm", draws = 30, seed = 1))
  s <- summary(r)
  expect_named(s, c("arm", "n", "mean", "lower", "upper"))
  expect_equal(s$arm, c("a", "b"))
  a <- r$mean_utility[r$arm == "a"]
  expect_equal(s$mean[1], mean(a))
  expect_equal(s$lower[1], unname(stats::quantile(a, 0.025)))
  expect_equal(s$upper[1], unname(stats::quantile(a, 0.975)))
})

test_that("as.matrix: rows x draws, groups x draws, error in mean mode", {
  nd <- res_data()
  ind <- quiet(predict_eq5d(nd, "UK", mode = "individual", draws = 7, seed = 1))
  m <- as.matrix(ind)
  expect_equal(dim(m), c(12L, 7L))
  expect_equal(m[3, 5], ind$utility[ind$.row == 3 & ind$.draw == 5])
  coh <- quiet(predict_eq5d(nd, "UK", mode = "cohort", group = "arm", draws = 7, seed = 1))
  mc <- as.matrix(coh)
  expect_equal(dim(mc), c(2L, 7L)); expect_equal(rownames(mc), c("a", "b"))
  expect_error(as.matrix(quiet(predict_eq5d(nd, "UK"))), "mean mode")
})
