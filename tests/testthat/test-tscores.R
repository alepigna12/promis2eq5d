fake_table <- function() do.call(rbind, lapply(seq_along(c("PF","ANX","DEP","FAT","SLP","SR","PI")), function(i)
  data.frame(domain = c("PF","ANX","DEP","FAT","SLP","SR","PI")[i], raw = 4:20, tscore = 10 * i + 4:20)))
test_that("raw sums map through the table", {
  d <- as.data.frame(as.list(setNames(rep(2, 28), ITEM_IDS)))
  ts <- p2e_tscores(d, fake_table())
  expect_equal(unlist(ts[1, ]), setNames(10 * 1:7 + 8, TSCORES))
  d$PFA11 <- NA; expect_true(is.na(p2e_tscores(d, fake_table())$T_PF))
  d$PFA11 <- 2; tb <- fake_table(); tb <- tb[!(tb$domain == "PF" & tb$raw == 8), ]
  expect_error(p2e_tscores(d, tb), "no table row.*PF")
  d$Sleep109 <- 5                                       # sleep = P5 + P6_1..3, so raw 5 + 2 + 2 + 2 = 11
  expect_equal(p2e_tscores(d, fake_table())$T_SLP, 10 * 5 + 11)
})
test_that(".standardise_input computes T-scores from a table", {
  d <- as.data.frame(as.list(setNames(rep(2, 28), ITEM_IDS))); d$Global07 <- 3; d$age <- 50; d$female <- 1
  expect_equal(.standardise_input(d, fake_table())$T_PF, 18)
})
test_that("computed T-scores replace supplied T_* columns with a warning", {
  d <- as.data.frame(as.list(setNames(rep(2, 28), ITEM_IDS))); d$Global07 <- 3; d$age <- 50; d$female <- 1
  d$T_PF <- 99
  expect_warning(s <- .standardise_input(d, fake_table()), "T_PF")
  expect_equal(s$T_PF, 18)
})
test_that("table and item columns are checked", {
  d <- as.data.frame(as.list(setNames(rep(2, 28), ITEM_IDS)))
  expect_error(p2e_tscores(d, fake_table()[, c("domain", "raw")]), "tscore")
  expect_error(p2e_tscores(d[, -1], fake_table()), "PFA11")
})
test_that("table validation: duplicates, unknown domains, NA, non-data-frame, no rows", {
  d <- as.data.frame(as.list(setNames(rep(2, 28), ITEM_IDS)))
  tb <- rbind(fake_table(), data.frame(domain = "PF", raw = 8, tscore = 99))
  expect_error(p2e_tscores(d, tb), "duplicate rows for PF raw 8")
  tb <- fake_table(); tb$domain[1] <- "XX"
  expect_error(p2e_tscores(d, tb), "XX.*PF")
  tb <- fake_table(); tb$tscore[3] <- NA
  expect_error(p2e_tscores(d, tb), "missing values.*tscore")
  expect_error(p2e_tscores(d, "x"), "data frame")
  expect_error(p2e_tscores(d, fake_table()[0, ]), "no rows")
})
test_that("a table with missing item columns is an error naming them", {
  d <- as.data.frame(as.list(setNames(rep(2, 28), ITEM_IDS))); d$Global07 <- 3; d$age <- 50; d$female <- 1
  d$PFA11 <- NULL
  expect_error(.standardise_input(d, fake_table()), "PFA11")
})
test_that("three rows keep their order across all seven T-scores", {
  d <- as.data.frame(as.list(setNames(rep(2, 28), ITEM_IDS)))[rep(1, 3), ]
  for (id in unname(ITEM_IDS)) d[[id]] <- c(2, 3, 4)
  ts <- p2e_tscores(d, fake_table())
  expect_equal(ts$T_PF, 10 * 1 + c(8, 12, 16))
  expect_equal(ts$T_SLP, 10 * 5 + c(8, 12, 16))
  expect_equal(ts$T_PI, 10 * 7 + c(8, 12, 16))
})
test_that("NA in one domain gives NA only for that domain", {
  d <- as.data.frame(as.list(setNames(rep(2, 28), ITEM_IDS)))[rep(1, 2), ]
  d[["EDANX01"]] <- c(NA, 2)
  ts <- p2e_tscores(d, fake_table())
  expect_true(is.na(ts$T_ANX[1])); expect_equal(ts$T_ANX[2], 10 * 2 + 8)
  expect_equal(ts$T_PF, c(18, 18)); expect_equal(ts$T_DEP, c(38, 38))
})
