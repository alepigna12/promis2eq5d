test_that("download verifies the checksum and installs into the cache", {
  src <- withr::local_tempfile(fileext = ".rds"); saveRDS(list(a = 1), src)
  withr::local_options(promis2eq5d.cache_dir = withr::local_tempdir(), promis2eq5d.allow_download = TRUE)
  p <- normalizePath(src, winslash = "/"); h <- .file_hash(src)
  local_mocked_bindings(.fit_manifest = function() data.frame(method = "R09_chain", set = "S3",
    url = paste0("file://", if (startsWith(p, "/")) "" else "/", p), sha256 = h[["sha256"]], md5 = h[["md5"]], bytes = file.size(src)))
  p2e_download("R09_chain", "S3", quiet = TRUE)
  expect_true(file.exists(file.path(.cache_dir(), "fits", "R09_chain_S3_any.rds")))
})
test_that("a bad checksum is an error and installs nothing", {
  src <- withr::local_tempfile(fileext = ".rds"); saveRDS(list(a = 1), src)
  withr::local_options(promis2eq5d.cache_dir = withr::local_tempdir(), promis2eq5d.allow_download = TRUE)
  p <- normalizePath(src, winslash = "/")
  local_mocked_bindings(.fit_manifest = function() data.frame(method = "R09_chain", set = "S3",
    url = paste0("file://", if (startsWith(p, "/")) "" else "/", p), sha256 = "0000", md5 = "0000", bytes = file.size(src)))
  expect_error(p2e_download("R09_chain", "S3", quiet = TRUE), "checksum")
  expect_false(file.exists(file.path(.cache_dir(), "fits", "R09_chain_S3_any.rds")))
})
test_that(".load_bundle errors non-interactively and offers the download interactively", {
  withr::local_options(promis2eq5d.fit_dir = tempfile(), promis2eq5d.cache_dir = withr::local_tempdir())
  expect_error(.load_bundle("R09_chain", "S3", NA), "p2e_download")
  # with mocks of interactive() = TRUE, utils::askYesNo = TRUE and p2e_download (records its call):
  # .load_bundle calls p2e_download once; with askYesNo = FALSE it errors with the hint.
})

# --- Remaining cases: mocked interactivity, consent and manifest ---

fake_bundle <- function() structure(list(method = "R09_chain"), class = "p2e_bundle")

test_that(".load_bundle calls p2e_download once when the user says yes, and loads the file", {
  cdir <- withr::local_tempdir()
  withr::local_options(promis2eq5d.fit_dir = tempfile(), promis2eq5d.cache_dir = cdir)
  local_mocked_bindings(.is_interactive = function() TRUE, .ask_yes_no = function(...) TRUE)
  n_calls <- 0L
  local_mocked_bindings(p2e_download = function(method = "R09_chain", set = "S3", quiet = FALSE) {
    n_calls <<- n_calls + 1L
    dir.create(file.path(cdir, "fits"), recursive = TRUE, showWarnings = FALSE)
    saveRDS(fake_bundle(), file.path(cdir, "fits", "R09_chain_S3_any.rds"))
    invisible(character())
  })
  b <- .load_bundle("R09_chain", "S3", NA)
  expect_equal(n_calls, 1L)
  expect_s3_class(b, "p2e_bundle")
})
test_that(".load_bundle errors with the hint when the user says no", {
  withr::local_options(promis2eq5d.fit_dir = tempfile(), promis2eq5d.cache_dir = withr::local_tempdir())
  local_mocked_bindings(.is_interactive = function() TRUE, .ask_yes_no = function(...) FALSE)
  n_calls <- 0L
  local_mocked_bindings(p2e_download = function(...) { n_calls <<- n_calls + 1L; invisible(character()) })
  expect_error(.load_bundle("R09_chain", "S3", NA), "p2e_download")
  expect_equal(n_calls, 0L)
})
test_that("no manifest row gives a clear 'not yet available' error", {
  withr::local_options(promis2eq5d.cache_dir = withr::local_tempdir(), promis2eq5d.allow_download = TRUE)
  local_mocked_bindings(.is_interactive = function() FALSE)
  expect_error(p2e_download("R09_chain", "S3", quiet = TRUE),
               "no download is available yet for R09_chain S3")
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
  expect_true(file.exists(file.path(cdir, "fits", "R09_chain_S3_any.rds")))
  expect_true(file.exists(file.path(cdir, "fits", "R09_chain_S1_any.rds")))
  expect_false(file.exists(file.path(cdir, "fits", "R09_chain_S4_any.rds")))
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
  expect_false(file.exists(file.path(cdir, "fits", "R09_chain_S3_any.rds")))
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
  f <- file.path(cdir, "fits", "R09_chain_S3_any.rds")
  expect_equal(readRDS(f), list(a = 2))
  expect_length(list.files(file.path(cdir, "fits"), all.files = TRUE, no.. = TRUE), 1L)
})
