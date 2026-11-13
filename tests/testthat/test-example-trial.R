test_that("p2e_example_trial structure", {
  d <- p2e_example_trial
  expect_identical(names(d), c("id", "arm", "month", unname(ITEM_IDS), PAIN_ID, "age", "female"))
  expect_equal(nrow(d), 900L)
  expect_true(all(table(d$id) == 3L)); expect_setequal(d$id, 1:300)
  expect_true(all(tapply(d$month, d$id, function(m) identical(sort(m), c(0L, 6L, 12L)))))
  expect_equal(as.vector(table(d$arm[d$month == 0])), c(150L, 150L))
  expect_setequal(d$arm, c("control", "treatment"))
  for (v in c("arm", "female")) expect_true(all(tapply(d[[v]], d$id, function(x) length(unique(x)) == 1L)), info = v)
  expect_true(all(tapply(d$age, d$id, function(x) all(abs(diff(x) - c(0.5, 0.5)) < 1e-9))))
})
test_that("p2e_example_trial values are in range", {
  d <- p2e_example_trial
  expect_false(anyNA(d))
  expect_true(all(vapply(d[unname(ITEM_IDS)], is.integer, logical(1))))
  expect_true(all(unlist(d[unname(ITEM_IDS)]) %in% 1:5))
  expect_true(all(d[[PAIN_ID]] %in% 0:10))
  expect_false(any(grepl("^T_|^E[1-5]$|country", names(d))))
})
test_that("p2e_example_trial predicts", {
  d <- p2e_example_trial[p2e_example_trial$month == 0, ]
  expect_message(r <- predict_eq5d(d, "UK"), "S1"); u <- r$utility
  expect_length(u, 300L); expect_true(all(is.finite(u))); expect_true(all(u <= 1))
})
test_that("the treatment effect shows", {
  d <- p2e_example_trial; pf <- unname(ITEM_IDS[1:4])
  chg <- function(a) {
    x <- d[d$arm == a, ]
    mean(unlist(x[x$month == 12, pf])) - mean(unlist(x[x$month == 0, pf]))
  }
  expect_gt(chg("treatment"), chg("control"))
})
