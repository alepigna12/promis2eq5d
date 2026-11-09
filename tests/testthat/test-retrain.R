retrain_data <- function() make_synthetic_data(400, seed = 1, countries = c("DE", "UK"))
retrain_item_data <- function() {
  d <- retrain_data()
  d[setdiff(names(d), TSCORES)]
}
quiet_retrain <- function(...) suppressMessages(.retrain_bundle(...))

test_that("a T-score set needs tscore_table", {
  withr::local_options(promis2eq5d.training_data = retrain_data(),
                       promis2eq5d.cache_dir = withr::local_tempdir())
  stand_in_bundle("D12_klapproth", "S3K", "UK")
  expect_error(.retrain_bundle("D12_klapproth", "S3K", p2e_value_set("USA")),
               "retrain on a T-score set needs tscore_table")
})

test_that("D03 retrain on an item set: target, boot, folds and caching", {
  cdir <- withr::local_tempdir()
  withr::local_options(promis2eq5d.training_data = retrain_item_data(), promis2eq5d.cache_dir = cdir)
  shipped <- stand_in_bundle("D03_clad", "S1", "UK")
  vs <- p2e_value_set("USA")
  rb <- quiet_retrain("D03", "S1", vs, boot = 2, mode = "individual", cache = FALSE)
  expect_s3_class(rb, "p2e_bundle")
  expect_equal(rb$target, "US"); expect_equal(rb$vs$code, "US")
  expect_length(rb$boot, 2); expect_length(rb$folds, 0)
  expect_identical(rb$fold_id, shipped$fold_id)
  expect_false(dir.exists(file.path(cdir, "retrain")))
  # individual draws work without fold fits
  nd <- retrain_item_data()[1:5, ]
  m <- .p2e_method("D03")
  expect_equal(dim(m$draw(rb$boot[[1]], nd, 3, key = 1:5, draw_type = NULL, lo = NULL)), c(5L, 3L))
  expect_equal(dim(m$draw(rb$point, nd, 3, key = 1:5, draw_type = NULL, lo = NULL)), c(5L, 3L))

  rbi <- quiet_retrain("D03", "S1", vs, boot = 0, interval = TRUE, mode = "individual", cache = FALSE)
  expect_length(rbi$folds, 3); expect_length(rbi$boot, 0)
  expect_length(rbi$oof$mu, 400)

  quiet_retrain("D03", "S1", vs, boot = 2, mode = "individual", cache = TRUE)
  expect_length(list.files(file.path(cdir, "retrain"), pattern = "[.]rds$"), 1)
})

test_that("a cached retrain is read back without refitting", {
  cdir <- withr::local_tempdir()
  withr::local_options(promis2eq5d.training_data = retrain_item_data(), promis2eq5d.cache_dir = cdir)
  stand_in_bundle("D03_clad", "S1", "UK")
  vs <- p2e_value_set("USA")
  msgs <- testthat::capture_messages(
    .retrain_bundle("D03", "S1", vs, boot = 1, mode = "individual", cache = TRUE))
  expect_true(any(grepl("point fit", msgs)))
  msgs2 <- testthat::capture_messages(
    r2 <- .retrain_bundle("D03", "S1", vs, boot = 1, mode = "individual", cache = TRUE))
  expect_true(any(grepl("Using cached refit", msgs2)))
  expect_false(any(grepl("point fit", msgs2)))
  expect_s3_class(r2, "p2e_bundle")
  # a mean-only bundle has a different key and is not reused for draws
  suppressMessages(.retrain_bundle("D03", "S1", vs, boot = 1, mode = "mean", cache = TRUE))
  expect_length(list.files(file.path(cdir, "retrain")), 2)
})

test_that("a mean-only D03 retrain can draw (it keeps its own residual pairs)", {
  withr::local_options(promis2eq5d.training_data = retrain_item_data(),
                       promis2eq5d.cache_dir = withr::local_tempdir())
  stand_in_bundle("D03_clad", "S1", "UK")
  rb <- quiet_retrain("D03", "S1", p2e_value_set("USA"), mode = "mean", cache = FALSE)
  expect_length(rb$boot, 0)
  nd <- retrain_item_data()[1:3, ]
  expect_false(is.null(rb$point$oof_resid))
  d <- .p2e_method("D03")$draw(rb$point, nd, 2, key = 1:3, draw_type = NULL, lo = NULL)
  expect_equal(dim(d), c(3L, 2L))
  expect_false(anyNA(d))
})

test_that("response methods and bad arguments error", {
  withr::local_options(promis2eq5d.training_data = retrain_item_data(),
                       promis2eq5d.cache_dir = withr::local_tempdir())
  vs <- p2e_value_set("USA")
  expect_error(.retrain_bundle("R06", "S1", vs), "response method")
  stand_in_bundle("D03_clad", "S1", "UK")
  expect_error(.retrain_bundle("D03", "S1", "US"), "p2e_vs")
  expect_error(.retrain_bundle("D03", "S1", vs, mode = "x"), "mode")
  expect_error(.retrain_bundle("D03", "S1", vs, interval = NA), "interval")
})

test_that(".training_data errors clearly when nothing is available", {
  withr::local_options(promis2eq5d.training_data = NULL)
  skip_if(exists("p2e_training", envir = asNamespace("promis2eq5d"), inherits = FALSE))
  expect_error(.training_data(), "training data")
})

test_that("fit_bundle: vs, fold_id and fit_folds", {
  train <- make_synthetic_data(200, 1)
  fid <- rep(1:3, length.out = 200)
  vs <- p2e_value_set("USA")
  b <- suppressMessages(.fit_bundle(toy_method(), train, "S4", "UK", list(), 3, 2, 1,
                                    vs = vs, fold_id = fid, fit_folds = FALSE, progress = TRUE))
  expect_equal(b$target, "US"); expect_identical(b$fold_id, fid)
  expect_length(b$folds, 0); expect_null(b$oof); expect_length(b$boot, 2)
  expect_error(.fit_bundle(toy_method(), train, "S4", "UK", list(), 3, 2, 1, fold_id = 1:5), "fold_id")
  msgs <- testthat::capture_messages(
    .fit_bundle(toy_method(), train, "S4", "UK", list(), 3, 1, 1, progress = TRUE))
  expect_equal(msgs, c("point fit\n", "fold 1/3\n", "fold 2/3\n", "fold 3/3\n", "boot 1/1\n"))
})

test_that("a T-score set retrains with a made-up table", {
  withr::local_options(promis2eq5d.training_data = retrain_data(),
                       promis2eq5d.cache_dir = withr::local_tempdir())
  stand_in_bundle("D12_klapproth", "S3K", "UK")
  doms <- c("PF", "ANX", "DEP", "FAT", "SLP", "SR", "PI")
  tb <- do.call(rbind, lapply(seq_along(doms), function(i)
    data.frame(domain = doms[i], raw = 0:40, tscore = 10 * i + (0:40) / 2)))
  rb <- quiet_retrain("D12_klapproth", "S3K", p2e_value_set("USA"), tscore_table = tb, cache = FALSE)
  expect_equal(rb$target, "US"); expect_length(rb$folds, 0); expect_length(rb$boot, 0)
})

test_that("a corrupt cache file warns, is refitted and overwritten", {
  cdir <- withr::local_tempdir()
  withr::local_options(promis2eq5d.training_data = retrain_item_data(), promis2eq5d.cache_dir = cdir)
  stand_in_bundle("D03_clad", "S1", "UK")
  vs <- p2e_value_set("USA")
  suppressMessages(.retrain_bundle("D03", "S1", vs, boot = 1, mode = "individual", cache = TRUE))
  path <- list.files(file.path(cdir, "retrain"), full.names = TRUE)
  expect_length(path, 1)
  saveRDS(list(a = 1), path)
  msgs <- testthat::capture_messages(
    expect_warning(r <- .retrain_bundle("D03", "S1", vs, boot = 1, mode = "individual", cache = TRUE),
                   "will be refitted"))
  expect_true(any(grepl("point fit", msgs)))
  expect_s3_class(r, "p2e_bundle")
  expect_s3_class(readRDS(path), "p2e_bundle")
})

test_that("T-score retrain replaces the input T_* columns from the table", {
  doms <- c("PF", "ANX", "DEP", "FAT", "SLP", "SR", "PI")
  tb <- do.call(rbind, lapply(seq_along(doms), function(i)
    data.frame(domain = doms[i], raw = 0:40, tscore = 10 * i + (0:40) / 2)))
  d <- retrain_data()
  dz <- d; dz[TSCORES] <- 0
  withr::local_options(promis2eq5d.cache_dir = withr::local_tempdir())
  withr::local_options(promis2eq5d.training_data = d)
  stand_in_bundle("D12_klapproth", "S3K", "UK")
  vs <- p2e_value_set("USA")
  a <- quiet_retrain("D12_klapproth", "S3K", vs, tscore_table = tb, cache = FALSE)
  withr::local_options(promis2eq5d.training_data = dz)
  z <- quiet_retrain("D12_klapproth", "S3K", vs, tscore_table = tb, cache = FALSE)
  nd <- d[1:20, ]; nd[names(.tscores_internal(nd, tb))] <- .tscores_internal(nd, tb)
  m <- .p2e_method("D12")
  expect_equal(m$point(z$point, nd), m$point(a$point, nd))
  expect_false(isTRUE(all.equal(m$point(a$point, nd), m$point(a$point, transform(nd, T_PF = 0, T_FAT = 0)))))
})

test_that("predict_eq5d reuses a saved refit without refitting", {
  cdir <- withr::local_tempdir()
  withr::local_options(promis2eq5d.training_data = retrain_item_data(), promis2eq5d.cache_dir = cdir)
  stand_in_bundle("D03_clad", "S1", "UK")
  suppressMessages(.retrain_bundle("D03", "S1", p2e_value_set("USA"), mode = "mean", cache = TRUE))
  expect_length(list.files(file.path(cdir, "retrain"), pattern = "[.]rds$"), 1)
  nd <- user_frame(make_synthetic_data(4, 2))
  msgs <- testthat::capture_messages(
    r <- predict_eq5d(nd, "USA", method = "D03", set = "S1", retrain = TRUE))
  expect_true(any(grepl("Using cached refit", msgs)))
  expect_false(any(grepl("point fit", msgs)))
  expect_length(r$utility, 4)
  expect_length(list.files(file.path(cdir, "retrain"), pattern = "[.]rds$"), 1)
})

test_that("T-score-only user data retrain a T-score set with a made-up table", {
  withr::local_options(promis2eq5d.training_data = retrain_data(),
                       promis2eq5d.cache_dir = withr::local_tempdir())
  stand_in_bundle("D12_klapproth", "S3K", "UK")
  doms <- c("PF", "ANX", "DEP", "FAT", "SLP", "SR", "PI")
  tb <- do.call(rbind, lapply(seq_along(doms), function(i)
    data.frame(domain = doms[i], raw = 0:40, tscore = 10 * i + (0:40) / 2)))
  nd <- make_synthetic_data(5, 2)[c(TSCORES, "AGE", "female")]
  names(nd)[names(nd) == "AGE"] <- "age"
  msgs <- testthat::capture_messages(
    r <- predict_eq5d(nd, "USA", method = "D12", retrain = TRUE, tscore_table = tb))
  expect_true(any(grepl("point fit", msgs)))
  expect_equal(attr(r, "p2e_info")$set, "S3K")
  expect_true(attr(r, "p2e_info")$retrain)
  expect_length(r$utility, 5); expect_false(anyNA(r$utility))
  # without the T-score columns the table cannot score incomplete items
  expect_error(suppressMessages(predict_eq5d(nd[c("age", "female")], "USA", method = "D12",
                                             retrain = TRUE, tscore_table = tb)),
               "neither all 28 items nor all T-score")
})
