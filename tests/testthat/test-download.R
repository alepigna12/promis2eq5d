test_that("download verifies the checksum and installs into the cache", {
  local_mocked_bindings(.is_interactive = function() FALSE)
  src <- withr::local_tempfile(fileext = ".rds"); saveRDS(list(a = 1), src)
  withr::local_options(promis2eq5d.cache_dir = withr::local_tempdir(), promis2eq5d.allow_download = TRUE)
  p <- normalizePath(src, winslash = "/"); h <- .file_hash(src)
  local_mocked_bindings(.fit_manifest = function() data.frame(method = "R09_chain", set = "S3",
    url = paste0("file://", if (startsWith(p, "/")) "" else "/", p), sha256 = h[["sha256"]], md5 = h[["md5"]], bytes = file.size(src)))
  p2e_download("R09_chain", "S3", quiet = TRUE)
  expect_true(file.exists(file.path(.cache_dir(), "fits", .bundle_file("R09_chain", "S3", NA))))
})
test_that("a bad checksum is an error and installs nothing", {
  local_mocked_bindings(.is_interactive = function() FALSE)
  src <- withr::local_tempfile(fileext = ".rds"); saveRDS(list(a = 1), src)
  withr::local_options(promis2eq5d.cache_dir = withr::local_tempdir(), promis2eq5d.allow_download = TRUE)
  p <- normalizePath(src, winslash = "/")
  local_mocked_bindings(.fit_manifest = function() data.frame(method = "R09_chain", set = "S3",
    url = paste0("file://", if (startsWith(p, "/")) "" else "/", p), sha256 = "0000", md5 = "0000", bytes = file.size(src)))
  expect_error(p2e_download("R09_chain", "S3", quiet = TRUE), "checksum")
  expect_false(file.exists(file.path(.cache_dir(), "fits", .bundle_file("R09_chain", "S3", NA))))
})
test_that(".load_bundle errors non-interactively and offers the download interactively", {
  local_mocked_bindings(.is_interactive = function() FALSE)
  withr::local_options(promis2eq5d.fit_dir = tempfile(), promis2eq5d.cache_dir = withr::local_tempdir())
  expect_error(.load_bundle("R09_chain", "S3", NA), "p2e_download")
})

# --- Remaining cases: mocked interactivity, consent and manifest ---

fake_bundle <- function() structure(list(method = "R09_chain", meta = list(format = 2L)), class = "p2e_bundle")

# A file:// URL for an existing local file
file_url <- function(path) {
  p <- normalizePath(path, winslash = "/")
  paste0("file://", if (startsWith(p, "/")) "" else "/", p)
}

one_row <- function(url, sha256, md5, bytes, set = "S3") {
  data.frame(method = "R09_chain", set = set, url = url, sha256 = sha256, md5 = md5,
             bytes = bytes, stringsAsFactors = FALSE)
}

test_that(".load_bundle asks once and calls p2e_download once when the user says yes", {
  cdir <- withr::local_tempdir()
  withr::local_options(promis2eq5d.fit_dir = tempfile(), promis2eq5d.cache_dir = cdir)
  local_mocked_bindings(.is_interactive = function() TRUE)
  local_mocked_bindings(.fit_manifest = function() one_row("file:///unused", NA, NA, 5e6))
  n_ask <- 0L
  local_mocked_bindings(.ask_yes_no = function(...) { n_ask <<- n_ask + 1L; TRUE })
  n_calls <- 0L
  consented <- NA
  local_mocked_bindings(p2e_download = function(method = "R09_chain", set = "S3", quiet = FALSE,
                                                .consented = FALSE) {
    n_calls <<- n_calls + 1L
    consented <<- .consented
    dir.create(file.path(cdir, "fits"), recursive = TRUE, showWarnings = FALSE)
    saveRDS(fake_bundle(), file.path(cdir, "fits", .bundle_file("R09_chain", "S3", NA)))
    invisible(character())
  })
  b <- .load_bundle("R09_chain", "S3", NA)
  expect_equal(n_ask, 1L)
  expect_equal(n_calls, 1L)
  expect_true(consented)
  expect_s3_class(b, "p2e_bundle")
})
test_that(".load_bundle asks once and errors with the hint when the user says no", {
  local_mocked_bindings(.is_interactive = function() TRUE)
  withr::local_options(promis2eq5d.fit_dir = tempfile(), promis2eq5d.cache_dir = withr::local_tempdir())
  local_mocked_bindings(.fit_manifest = function() one_row("file:///unused", NA, NA, 5e6))
  n_ask <- 0L
  local_mocked_bindings(.ask_yes_no = function(...) { n_ask <<- n_ask + 1L; FALSE })
  n_calls <- 0L
  local_mocked_bindings(p2e_download = function(...) { n_calls <<- n_calls + 1L; invisible(character()) })
  expect_error(.load_bundle("R09_chain", "S3", NA), "p2e_download")
  expect_equal(n_ask, 1L)
  expect_equal(n_calls, 0L)
})
test_that(".load_bundle without a manifest row errors with the hint and asks nothing", {
  local_mocked_bindings(.is_interactive = function() TRUE)
  withr::local_options(promis2eq5d.fit_dir = tempfile(), promis2eq5d.cache_dir = withr::local_tempdir())
  local_mocked_bindings(.fit_manifest = function() data.frame(method = character(), set = character(),
    url = character(), sha256 = character(), md5 = character(), bytes = numeric()))
  n_ask <- 0L
  local_mocked_bindings(.ask_yes_no = function(...) { n_ask <<- n_ask + 1L; TRUE })
  expect_error(.load_bundle("R09_chain", "S3", NA),
               "no download is available yet for R09_chain S3.*p2e_download\\('R09_chain', set = 'S3'\\)")
  expect_equal(n_ask, 0L)
})
test_that("no manifest row gives a clear 'not yet available' error", {
  local_mocked_bindings(.is_interactive = function() FALSE)
  withr::local_options(promis2eq5d.cache_dir = withr::local_tempdir(), promis2eq5d.allow_download = TRUE)
  local_mocked_bindings(.fit_manifest = function() data.frame(method = character(), set = character(),
    url = character(), sha256 = character(), md5 = character(), bytes = numeric()))
  expect_error(p2e_download("R09_chain", "S3", quiet = TRUE),
               "no download is available yet for R09_chain S3")
})
test_that("duplicate manifest rows are a distinct error", {
  local_mocked_bindings(.is_interactive = function() FALSE)
  withr::local_options(promis2eq5d.cache_dir = withr::local_tempdir(), promis2eq5d.allow_download = TRUE)
  local_mocked_bindings(.fit_manifest = function() rbind(one_row("file:///a", NA, "x", 1),
                                                        one_row("file:///b", NA, "y", 1)))
  expect_error(p2e_download("R09_chain", "S3", quiet = TRUE),
               "the manifest has duplicate rows for R09_chain S3")
})
test_that("set = 'all' downloads every set that has a manifest row", {
  src <- withr::local_tempfile(fileext = ".rds"); saveRDS(list(a = 1), src)
  cdir <- withr::local_tempdir()
  withr::local_options(promis2eq5d.cache_dir = cdir, promis2eq5d.allow_download = TRUE)
  local_mocked_bindings(.is_interactive = function() FALSE)
  p <- normalizePath(src, winslash = "/"); h <- .file_hash(src)
  u <- paste0("file://", if (startsWith(p, "/")) "" else "/", p)
  local_mocked_bindings(.fit_manifest = function() data.frame(method = "R09_chain", set = c("S3", "S1"),
    url = c(u, u), sha256 = h[["sha256"]], md5 = h[["md5"]], bytes = c(file.size(src), file.size(src))))
  p2e_download("R09_chain", "all", quiet = TRUE)
  expect_true(file.exists(file.path(cdir, "fits", .bundle_file("R09_chain", "S3", NA))))
  expect_true(file.exists(file.path(cdir, "fits", .bundle_file("R09_chain", "S1", NA))))
  expect_false(file.exists(file.path(cdir, "fits", .bundle_file("R09_chain", "S4", NA))))
})
test_that("a non-interactive download without the option errors with the command and installs nothing", {
  src <- withr::local_tempfile(fileext = ".rds"); saveRDS(list(a = 1), src)
  cdir <- withr::local_tempdir()
  withr::local_options(promis2eq5d.cache_dir = cdir, promis2eq5d.allow_download = NULL)
  local_mocked_bindings(.is_interactive = function() FALSE)
  p <- normalizePath(src, winslash = "/"); h <- .file_hash(src)
  local_mocked_bindings(.fit_manifest = function() data.frame(method = "R09_chain", set = "S3",
    url = paste0("file://", if (startsWith(p, "/")) "" else "/", p), sha256 = h[["sha256"]], md5 = h[["md5"]], bytes = file.size(src)))
  expect_error(p2e_download("R09_chain", "S3", quiet = TRUE),
               "options\\(promis2eq5d.allow_download = TRUE\\)", fixed = FALSE)
  expect_false(file.exists(file.path(cdir, "fits", .bundle_file("R09_chain", "S3", NA))))
})
test_that("a repeated download replaces the installed file", {
  src <- withr::local_tempfile(fileext = ".rds"); saveRDS(list(a = 2), src)
  cdir <- withr::local_tempdir()
  withr::local_options(promis2eq5d.cache_dir = cdir, promis2eq5d.allow_download = TRUE)
  local_mocked_bindings(.is_interactive = function() FALSE)
  p <- normalizePath(src, winslash = "/"); h <- .file_hash(src)
  local_mocked_bindings(.fit_manifest = function() data.frame(method = "R09_chain", set = "S3",
    url = paste0("file://", if (startsWith(p, "/")) "" else "/", p), sha256 = h[["sha256"]], md5 = h[["md5"]], bytes = file.size(src)))
  p2e_download("R09_chain", "S3", quiet = TRUE)
  p2e_download("R09_chain", "S3", quiet = TRUE)
  f <- file.path(cdir, "fits", .bundle_file("R09_chain", "S3", NA))
  expect_equal(readRDS(f), list(a = 2))
  expect_length(list.files(file.path(cdir, "fits"), all.files = TRUE, no.. = TRUE), 1L)
})
test_that("the download timeout is raised during the transfer and restored afterwards", {
  src <- withr::local_tempfile(fileext = ".rds"); saveRDS(list(a = 1), src)
  cdir <- withr::local_tempdir()
  withr::local_options(promis2eq5d.cache_dir = cdir, promis2eq5d.allow_download = TRUE, timeout = 60)
  local_mocked_bindings(.is_interactive = function() FALSE)
  h <- .file_hash(src)
  local_mocked_bindings(.fit_manifest = function() one_row(file_url(src), h[["sha256"]], h[["md5"]],
                                                           file.size(src)))
  seen <- NA_real_
  local_mocked_bindings(.download_file = function(url, dest, quiet) {
    seen <<- getOption("timeout")
    file.copy(src, dest, overwrite = TRUE)
    invisible(0L)
  })
  p2e_download("R09_chain", "S3", quiet = TRUE)
  expect_gte(seen, 3600)
  expect_equal(getOption("timeout"), 60)
  expect_true(file.exists(file.path(cdir, "fits", .bundle_file("R09_chain", "S3", NA))))
})
test_that("a failed download errors with context, installs nothing and leaves no partial files", {
  cdir <- withr::local_tempdir()
  withr::local_options(promis2eq5d.cache_dir = cdir, promis2eq5d.allow_download = TRUE)
  local_mocked_bindings(.is_interactive = function() FALSE)
  d <- withr::local_tempdir()
  missing_url <- paste0("file://", normalizePath(d, winslash = "/"), "/does-not-exist.rds")
  local_mocked_bindings(.fit_manifest = function() one_row(missing_url, "0000", "0000", 1))
  expect_error(p2e_download("R09_chain", "S3", quiet = TRUE), "download of R09_chain S3 failed")
  expect_false(file.exists(file.path(cdir, "fits", .bundle_file("R09_chain", "S3", NA))))
  expect_length(list.files(cdir, recursive = TRUE, all.files = TRUE, pattern = "\\.part$"), 0L)
})
