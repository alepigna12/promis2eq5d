# The fits shipped in inst/fits (R06, D03 and D12; R09 is downloaded)

test_that("every shipped bundle loads, has format 2 and the documented sizes", {
  files <- list.files(system.file("fits", package = "promis2eq5d"), full.names = TRUE)
  expect_length(files, 19L)
  digest <- .data_digest(.training_data())
  for (f in files) {
    b <- readRDS(f)
    expect_s3_class(b, "p2e_bundle")
    expect_identical(b$meta$format, 2L)
    expect_equal(b$meta$K, 10L)
    expect_equal(b$meta$n, 4512L)
    expect_identical(b$meta$train_digest, digest)
    expect_true(b$meta$B %in% c(50L, 200L))
    expect_false(is.na(b$meta$repo_commit))
    expect_true(length(b$boot) == b$meta$B || b$boot_attempted == 2L * b$meta$B)
  }
})

test_that("the defaults predict on training rows in every mode", {
  withr::local_options(promis2eq5d.fit_dir = NULL)
  d <- as.data.frame(.training_data()[1:5, ])
  cols <- p2e_columns()
  cols <- cols[cols$type == "item" | cols$column %in% c("age", "female", "Global07"), ]
  pro <- as.data.frame(lapply(seq_len(nrow(cols)), function(i) d[[cols$internal[i]]]))
  names(pro) <- cols$column
  lo <- .vs_min(p2e_value_set("UK"))
  ok <- function(x) all(is.finite(x) & x >= lo - 1e-8 & x <= 1 + 1e-8)
  for (fam in c("response", "direct")) {
    m <- suppressMessages(predict_eq5d(pro, "UK", family = fam))
    expect_equal(nrow(m), 5L)
    expect_true(ok(m$utility))
    ind <- suppressMessages(predict_eq5d(pro, "UK", family = fam, mode = "individual",
                                         draws = 20, seed = 1))
    expect_equal(nrow(ind), 5L * 20L)
    expect_true(ok(ind$utility))
    coh <- suppressMessages(predict_eq5d(pro, "UK", family = fam, mode = "cohort",
                                         draws = 20, seed = 1))
    expect_true(ok(as.matrix(coh)))
  }
})
