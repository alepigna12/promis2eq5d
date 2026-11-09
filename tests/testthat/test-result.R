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

test_that("print, summary and as.matrix work on a column subset", {
  nd <- res_data(); nd$pid <- paste0("p", 1:12)
  r <- quiet(predict_eq5d(nd, "UK", id = "pid"))
  sub <- r[c("pid", "utility")]
  expect_s3_class(sub, "p2e_result")
  out <- capture.output(print(sub))
  expect_true(any(grepl("utility", out)))
  expect_false(any(grepl("method:", out)))
  expect_s3_class(summary(sub), "table")            # the plain data-frame summary
  expect_equal(dim(as.matrix(sub)), c(12L, 2L))
})

test_that("summary and as.matrix use the recorded id and group columns", {
  nd <- res_data(); nd$pid <- paste0("p", 1:12)
  ind <- quiet(predict_eq5d(nd, "UK", mode = "individual", draws = 5, seed = 1, id = "pid"))
  coh <- quiet(predict_eq5d(nd, "UK", mode = "cohort", group = "arm", draws = 5, seed = 1))
  expect_equal(attr(ind, "p2e_cols"), list(id = "pid", group = NULL))
  expect_equal(attr(coh, "p2e_cols"), list(id = NULL, group = "arm"))
  s_ind <- summary(ind); s_coh <- summary(coh); m_coh <- as.matrix(coh)
  ind$extra <- 1; coh$extra <- 1
  expect_equal(summary(ind), s_ind)
  expect_equal(summary(coh), s_coh)
  expect_equal(as.matrix(coh), m_coh)
})

test_that("draw modes report the bootstrap counts in print and summary", {
  r <- quiet(predict_eq5d(res_data(), "UK", mode = "individual", draws = 4, seed = 1))
  i <- attr(r, "p2e_info")
  expect_equal(i$boot_attempted, 4L); expect_equal(i$boot_failed, 0L)
  out <- capture.output(print(r))
  expect_true(any(grepl("bootstrap fits: 4 of 4 converged", out, fixed = TRUE)))
  expect_equal(attr(summary(r), "boot_attempted"), 4L)
  expect_equal(attr(summary(r), "boot_failed"), 0L)
  # failures are subtracted
  attr(r, "p2e_info")$boot_failed <- 1L
  expect_true(any(grepl("bootstrap fits: 3 of 4 converged", capture.output(print(r)), fixed = TRUE)))
  co <- quiet(predict_eq5d(res_data(), "UK", mode = "cohort", group = "arm", draws = 4, seed = 1))
  expect_equal(attr(summary(co), "boot_attempted"), 4L)
  # mean mode: no line, no attributes
  rm <- quiet(predict_eq5d(res_data(), "UK"))
  expect_false(any(grepl("bootstrap fits", capture.output(print(rm)))))
  expect_null(attr(summary(rm), "boot_attempted"))
})
