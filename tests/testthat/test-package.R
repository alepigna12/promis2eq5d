test_that("package metadata is as specified", {
  d <- utils::packageDescription("promis2eq5d")
  expect_equal(d$License, "MIT + file LICENSE")
  expect_match(d$Title, "'PROMIS'", fixed = TRUE)
  expect_match(d$Title, "'EQ-5D'", fixed = TRUE)
  expect_match(d$Imports, "eq5d \\(>= 0\\.17\\.0\\)")
})

test_that("Description cites the training data DOI", {
  d <- utils::packageDescription("promis2eq5d")
  expect_true(nzchar(d$Description))
  expect_match(d$Description, "doi:10.7910/DVN/FNYYQR", fixed = TRUE)
})
