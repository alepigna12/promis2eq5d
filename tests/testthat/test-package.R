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

test_that("citation() lists the package, the R09 fits and the data", {
  ct <- citation("promis2eq5d")
  expect_length(ct, 3)
  txt <- paste(format(ct, style = "text"), collapse = " ")
  expect_match(txt, "10.5281/zenodo.24199380", fixed = TRUE)
  expect_match(txt, "10.7910/DVN/FNYYQR", fixed = TRUE)
})
