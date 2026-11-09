use_stand_in_fits()
quiet <- function(expr) suppressMessages(expr)
# errors raised after the set message: keep the output clean
expect_error <- function(object, ...) testthat::expect_error(suppressMessages(object), ...)
# the messages of a call must include one matching regex; the others are muffled
expect_msg <- function(expr, regexp) {
  msgs <- testthat::capture_messages(expr)
  testthat::expect_match(msgs, regexp, all = FALSE)
}

# Collect the warnings of an expression without printing them
warnings_of <- function(expr) {
  w <- character()
  withCallingHandlers(expr, warning = function(c) {
    w <<- c(w, conditionMessage(c)); invokeRestart("muffleWarning")
  })
  w
}

test_that("family/method resolution covers 0, 1 and 2 arguments", {
  skip_if_not_installed("xgboost")                 # method = "chain" is R09
  nd <- user_frame(make_synthetic_data(5, 2))
  expect_equal(attr(quiet(predict_eq5d(nd, "UK")), "p2e_info")$method, "R06_mvop_dim")
  expect_equal(attr(quiet(predict_eq5d(nd, "UK", family = "direct")), "p2e_info")$method, "D03_clad")
  expect_equal(attr(quiet(predict_eq5d(nd, "UK", method = "chain")), "p2e_info")$family, "response")
  expect_error(predict_eq5d(nd, "UK", family = "direct", method = "R06"), "R06_mvop_dim is a response method")
  expect_error(predict_eq5d(nd, "UK", family = "both"), "family")
})

test_that("cohort draws are the per-draw means of individual draws", {
  nd <- user_frame(make_synthetic_data(12, 2)); nd$arm <- rep(c("b", "a"), 6)   # first group is not first alphabetically
  ind <- quiet(predict_eq5d(nd, "UK", mode = "individual", draws = 8, seed = 1))
  coh <- quiet(predict_eq5d(nd, "UK", mode = "cohort", group = "arm", draws = 8, seed = 1))
  manual <- aggregate(utility ~ .draw, data = transform(ind, arm = nd$arm[.row])[nd$arm[ind$.row] == "a", ], mean)
  expect_equal(coh$mean_utility[coh$arm == "a"], manual$utility)
  expect_named(coh, c("arm", ".draw", "n", "mean_utility"))
  expect_equal(unique(coh$n[coh$arm == "a"]), 6L)
  expect_equal(unique(coh$arm), c("b", "a"))   # ordered by first appearance
})

test_that("crosswalk and response-retrain routes are announced", {
  nd <- user_frame(make_synthetic_data(5, 2))
  expect_msg(predict_eq5d(nd, "USA", family = "direct", draws = 4), "iso crosswalk from (FR|DE|UK) to USA")
  expect_msg(predict_eq5d(nd, "UK", retrain = TRUE), "no effect for response methods")
  expect_error(predict_eq5d(nd, "UK", family = "direct", dimensions = TRUE), "need a response method")
  expect_error(predict_eq5d(nd, "UK", mode = "individual", family = "direct", dimensions = TRUE),
               "need a response method")
  expect_msg(predict_eq5d(nd, "UK"), "Using predictor set S1 \\(")
  expect_error(predict_eq5d(nd, "USA", family = "direct", crosswalk_from = "JP"), "FR, DE, UK")
  expect_msg(predict_eq5d(nd, "USA", family = "direct", crosswalk_from = "FR"), "crosswalk from FR to USA")
})

test_that("tibbles, grouped tibbles, one row and NA ids", {
  skip_if_not_installed("tibble")
  nd <- user_frame(make_synthetic_data(6, 2)); nd$pid <- rep(1:3, 2)
  expect_equal(quiet(predict_eq5d(tibble::as_tibble(nd), "UK"))$utility, quiet(predict_eq5d(nd, "UK"))$utility)
  expect_s3_class(quiet(predict_eq5d(nd[1, ], "UK", mode = "cohort", draws = 3)), "p2e_result")
  nd$pid[2] <- NA
  expect_error(predict_eq5d(nd, "UK", mode = "individual", id = "pid", within_id = "common"), "pid.*missing")
})

test_that("a grouped tibble gives a plain data-frame result", {
  skip_if_not_installed("tibble")
  nd <- user_frame(make_synthetic_data(6, 2))
  g <- tibble::as_tibble(nd); class(g) <- c("grouped_df", "tbl_df", "tbl", "data.frame")
  r <- quiet(predict_eq5d(g, "UK"))
  expect_identical(class(r), c("p2e_result", "data.frame"))
  expect_equal(r$utility, quiet(predict_eq5d(nd, "UK"))$utility)
})

test_that("D03 default statistic is the mean, and differs from the median", {
  nd <- user_frame(make_synthetic_data(5, 2))
  mn <- quiet(predict_eq5d(nd, "UK", family = "direct"))$utility
  md <- quiet(predict_eq5d(nd, "UK", family = "direct", statistic = "median"))$utility
  expect_false(isTRUE(all.equal(mn, md)))
  expect_equal(mn, quiet(predict_eq5d(nd, "UK", family = "direct"))$utility)  # fixed MC seed
})

test_that("mean mode columns, interval and dimensions", {
  nd <- user_frame(make_synthetic_data(8, 2)); nd$pid <- letters[1:8]
  r <- quiet(predict_eq5d(nd, "UK", interval = TRUE, dimensions = TRUE, id = "pid"))
  expect_named(r, c(".row", "pid", "utility", "lower", "upper", "MO", "SC", "UA", "PD", "AD"))
  expect_equal(r$.row, 1:8)
  expect_true(all(r$lower <= r$utility & r$utility <= r$upper))
  expect_true(all(as.matrix(r[c("MO", "SC", "UA", "PD", "AD")]) %in% 1:5))
  # a direct method (D12) gets an interval from the shipped fold fits
  d <- quiet(predict_eq5d(nd, "UK", method = "D12", interval = TRUE))
  expect_true(all(d$lower <= d$utility & d$utility <= d$upper))
  expect_identical(class(r), c("p2e_result", "data.frame"))
})

test_that("individual mode: long format, dimensions, one row", {
  nd <- user_frame(make_synthetic_data(5, 2))
  r <- quiet(predict_eq5d(nd, "UK", mode = "individual", draws = 6, dimensions = TRUE, seed = 3))
  expect_equal(nrow(r), 30)
  expect_named(r, c(".row", ".draw", "utility", "MO", "SC", "UA", "PD", "AD"))
  expect_equal(r$utility[1:5], .score(as.matrix(r[1:5, c("MO", "SC", "UA", "PD", "AD")]), p2e_value_set("UK")))
  one <- quiet(predict_eq5d(nd[1, ], "UK", mode = "individual", draws = 4, seed = 3))
  expect_equal(nrow(one), 4)
  m1 <- quiet(predict_eq5d(nd[1, ], "UK"))
  expect_equal(nrow(m1), 1)
  d1 <- quiet(predict_eq5d(nd[1, ], "UK", family = "direct"))
  expect_equal(nrow(d1), 1)
})

test_that("seed makes draws reproducible and within_id = common correlates rows", {
  nd <- user_frame(make_synthetic_data(6, 2)); nd$pid <- rep(1:3, 2)
  a <- quiet(predict_eq5d(nd, "UK", mode = "individual", draws = 5, seed = 9))
  b <- quiet(predict_eq5d(nd, "UK", mode = "individual", draws = 5, seed = 9))
  expect_equal(a$utility, b$utility)
  # identical rows with a common id give identical draws
  two <- nd[c(1, 1), ]; two$pid <- 1
  ind <- quiet(predict_eq5d(two, "UK", mode = "individual", draws = 20, seed = 2, id = "pid"))
  com <- quiet(predict_eq5d(two, "UK", mode = "individual", draws = 20, seed = 2, id = "pid", within_id = "common"))
  expect_false(isTRUE(all.equal(ind$utility[ind$.row == 1], ind$utility[ind$.row == 2])))
  expect_equal(com$utility[com$.row == 1], com$utility[com$.row == 2])
  expect_equal(attr(com, "p2e_info")$within_id, "common")
})

test_that("group NA, missing group or id columns, and unused factor levels", {
  nd <- user_frame(make_synthetic_data(6, 2)); nd$arm <- rep(c("a", "b"), 3)
  g <- nd; g$arm[3] <- NA
  expect_error(predict_eq5d(g, "UK", mode = "cohort", group = "arm", draws = 3), "arm.*missing")
  expect_error(predict_eq5d(nd, "UK", mode = "cohort", group = "nope", draws = 3), "nope")
  expect_error(predict_eq5d(nd, "UK", id = "nope"), "nope")
  nd$f <- factor(nd$arm, levels = c("a", "b", "unused"))
  r <- quiet(predict_eq5d(nd, "UK", mode = "cohort", group = "f", draws = 3))
  expect_equal(nrow(r), 6)
  expect_false("unused" %in% r$f)
  expect_equal(levels(r$f), c("a", "b"))
  all <- quiet(predict_eq5d(nd, "UK", mode = "cohort", draws = 3))
  expect_equal(unique(all$n), 6L)
  two <- quiet(predict_eq5d(transform(nd, t = rep(1:2, each = 3)), "UK", mode = "cohort",
                            group = c("arm", "t"), draws = 2))
  expect_named(two, c("arm", "t", ".draw", "n", "mean_utility"))
})

test_that("as_published warns exactly once per call", {
  nd <- user_frame(make_synthetic_data(5, 2))
  w <- warnings_of(quiet(predict_eq5d(nd, "UK", method = "D12", draw_type = "as_published")))
  expect_length(w, 1); expect_match(w, "as_published")
  w <- warnings_of(quiet(predict_eq5d(nd, "USA", method = "D03", draw_type = "as_published", draws = 4)))
  expect_length(w, 1)
  w <- warnings_of(quiet(predict_eq5d(nd, "UK", method = "D03", draw_type = "as_published",
                                      mode = "individual", draws = 4)))
  expect_length(w, 1)
  w <- warnings_of(quiet(predict_eq5d(nd, "UK", method = "D12", draw_type = "as_published",
                                      mode = "cohort", draws = 4)))
  expect_length(w, 1)
  expect_length(warnings_of(quiet(predict_eq5d(nd, "UK", method = "D12"))), 0)
})

test_that("result attributes record the call", {
  nd <- user_frame(make_synthetic_data(5, 2))
  r <- quiet(predict_eq5d(nd, "USA", family = "direct", mode = "individual", draws = 4, seed = 5))
  i <- attr(r, "p2e_info")
  expect_equal(i$method, "D03_clad"); expect_equal(i$family, "direct")
  expect_equal(i$set, "S3"); expect_equal(i$mode, "individual"); expect_equal(i$draws, 4)
  expect_equal(i$value_set$code, "US"); expect_true(i$crosswalk %in% c("FR", "DE", "UK"))
  expect_false(i$retrain); expect_equal(i$seed, 5)
  expect_equal(i$version, utils::packageVersion("promis2eq5d"))
  expect_equal(i$B, 4)
})

test_that("crosswalked direct mean is the mean of the crosswalked draws", {
  nd <- user_frame(make_synthetic_data(5, 2))
  usa <- p2e_value_set("USA")
  r <- quiet(predict_eq5d(nd, "USA", family = "direct", crosswalk_from = "UK"))
  expect_true(all(r$utility <= 1 & r$utility >= min(usa$u)))
  # the same Monte Carlo mean by hand: exact, so the default tolerance applies
  b <- .load_bundle("D03_clad", "S3", "UK")
  knots <- .fit_crosswalk(b$vs, usa, b$state_freq)
  std <- .standardise_input(nd)
  manual <- .with_mc_seed(rowMeans(.draw_matrix(b, std, usa, N_MC_MEAN, 1:5, "corrected",
                                                 crosswalk = knots, fits = list(b$point))$u))
  expect_equal(r$utility, manual)
})

test_that("the D03 Monte Carlo mean does not depend on the RNG kind", {
  nd <- user_frame(make_synthetic_data(5, 2))
  base_uk <- quiet(predict_eq5d(nd, "UK", family = "direct"))$utility
  base_us <- quiet(predict_eq5d(nd, "USA", family = "direct", crosswalk_from = "UK"))$utility
  old <- RNGkind()
  withr::defer(do.call(RNGkind, as.list(old)))
  RNGkind("L'Ecuyer-CMRG")
  expect_identical(quiet(predict_eq5d(nd, "UK", family = "direct"))$utility, base_uk)
  expect_identical(quiet(predict_eq5d(nd, "USA", family = "direct", crosswalk_from = "UK"))$utility, base_us)
  expect_identical(RNGkind()[1], "L'Ecuyer-CMRG")   # the caller's kind is left alone
})

test_that("the D03 Monte Carlo mean is computed in row blocks", {
  nd <- user_frame(make_synthetic_data(5, 2))
  whole <- quiet(predict_eq5d(nd, "UK", family = "direct"))$utility
  local_mocked_bindings(MC_MEAN_BLOCK = 100L)            # n <= block: identical
  expect_identical(quiet(predict_eq5d(nd, "UK", family = "direct"))$utility, whole)
  local_mocked_bindings(MC_MEAN_BLOCK = 2L)
  a <- quiet(predict_eq5d(nd, "UK", family = "direct"))$utility
  expect_length(a, 5)
  expect_identical(quiet(predict_eq5d(nd, "UK", family = "direct"))$utility, a)  # deterministic
  expect_true(all(a >= 0 & a <= 1))
})

test_that("ignored arguments are announced and level is checked up front", {
  nd <- user_frame(make_synthetic_data(5, 2)); nd$arm <- rep(c("a", "b"), length.out = 5)
  expect_msg(predict_eq5d(nd, "UK", family = "direct", retrain = TRUE),
             "retrain has no effect: UK has a shipped fit")
  w <- warnings_of(quiet(predict_eq5d(nd, "UK", group = "arm")))
  expect_length(w, 1); expect_match(w, "group is only used in cohort mode")
  expect_error(predict_eq5d(nd, "UK", interval = TRUE, level = 1.5), "level")
  expect_error(predict_eq5d(nd, "UK", interval = TRUE, level = c(.8, .9)), "level")
})

test_that("a missing predict package is an error with an install hint", {
  skip_if_not_installed("xgboost")
  nd <- user_frame(make_synthetic_data(5, 2))
  local_mocked_bindings(requireNamespace = function(package, ..., quietly = FALSE) FALSE, .package = "base")
  expect_error(predict_eq5d(nd, "UK", method = "R09"), "xgboost.*install")
})

test_that("bad arguments are errors", {
  nd <- user_frame(make_synthetic_data(5, 2))
  expect_error(predict_eq5d("x", "UK"), "data frame")
  expect_error(predict_eq5d(nd[0, ], "UK"), "no rows")
  expect_error(predict_eq5d(nd, "UK", draws = 0, mode = "individual"), "draws")
  expect_error(predict_eq5d(nd, "UK", mode = "bad"), "should be one of")
  expect_error(predict_eq5d(nd, "UK", seed = "a"), "seed")
})
