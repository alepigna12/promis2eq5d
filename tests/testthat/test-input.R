test_that("item IDs map to internal names one to one", {
  expect_length(ITEM_IDS, 28); expect_setequal(names(ITEM_IDS), ITEMS)
  expect_equal(unname(ITEM_IDS[c("P1_1", "P5", "P7_4", "P8_4")]), c("PFA11", "Sleep109", "SRPPER46_CaPS", "PAININ34"))
})
test_that("standardise maps names, converts sex, keeps row order", {
  d <- data.frame(T_PF = 40, T_ANX = 55, T_DEP = 50, T_FAT = 52, T_SLP = 48, T_SR = 45, T_PI = 60,
                  Global07 = 4, age = "52", sex = "F")
  expect_error(.standardise_input(d), "age.*numeric")
  d$age <- 52
  s <- .standardise_input(d)
  expect_equal(s$P9, 4); expect_equal(s$AGE, 52); expect_equal(s$female, 1)
  expect_setequal(.available_sets(s), c("S3", "S4", "S3K"))
})
test_that("invalid values are errors naming the column", {
  d <- as.data.frame(as.list(setNames(rep(3, 28), ITEM_IDS))); d$Global07 <- 11
  expect_error(.standardise_input(d), "Global07.*0-10")
  d$Global07 <- 3; d$PFA11 <- factor("3")
  expect_error(.standardise_input(d), "PFA11.*factor")
  d$PFA11 <- 6; expect_error(.standardise_input(d), "PFA11.*1-5")
  d$PFA11 <- NA; expect_error(.standardise_input(d), "PFA11.*missing.*row 1")
})
test_that("set choice follows the preferred order and honours overrides", {
  items <- as.data.frame(as.list(setNames(rep(3, 28), ITEM_IDS))); items$Global07 <- 2
  s <- .standardise_input(items)
  expect_equal(.choose_set("R09_chain", s), "S2")           # S3, S1, S4 unavailable
  expect_equal(.choose_set("R06_mvop_dim", s), "S2")
  expect_error(.choose_set("D12_klapproth", s), "T-scores.*age.*sex")
  expect_error(.choose_set("R06_mvop_dim", s, set = "S3"), "S3.*not supported by the data")
})
test_that("columns table, age warning and aliases", {
  pc <- p2e_columns(); expect_setequal(pc$column[pc$type == "item"], unname(ITEM_IDS))
  d <- data.frame(T_PF = 40, T_ANX = 55, T_DEP = 50, T_FAT = 52, T_SLP = 48, T_SR = 45, T_PI = 60, P9 = 4, AGE = 15, female = 0)
  expect_warning(s <- .standardise_input(d), "18-100")
  expect_equal(c(s$P9, s$AGE, s$female), c(4, 15, 0))
  d$female <- NULL; d$sex <- "male"; expect_equal(suppressWarnings(.standardise_input(d))$female, 0)
})
