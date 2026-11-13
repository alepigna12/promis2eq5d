#' Download the large fit bundles
#'
#' Downloads the fitted model bundles that `predict_eq5d()` needs for the
#' R09_chain method, verifies each file against the published checksum and
#' installs it in the cache, from where it is found automatically. The files
#' are hosted on Zenodo, \doi{10.5281/zenodo.24199380}.
#'
#' @section Consent:
#' In an interactive session the function asks before downloading and states
#' the total size in MB. In a non-interactive session it downloads only when
#' `options(promis2eq5d.allow_download = TRUE)` is set; otherwise it stops
#' with the command to run.
#'
#' @section Cache location:
#' Files are installed in `file.path(tools::R_user_dir("promis2eq5d", "cache"),
#' "fits")`, unless `options(promis2eq5d.cache_dir)` is set. Each file is
#' named `<method>_<set>_any.rds`.
#'
#' @section Checksum:
#' Each file is checked against its SHA-256 checksum (computed with
#' `tools::sha256sum()`, available from R 4.5.0). On older versions of R the
#' MD5 checksum is used instead. A mismatch is an error and nothing is
#' installed.
#'
#' @section Timeout:
#' The download timeout is raised to at least one hour while a file is
#' transferred, because the fits are large, and restored afterwards.
#'
#' @param method The response method to download. Only `"R09_chain"` has
#'   downloadable fits.
#' @param set The predictor set to download: `"S3"`, `"S1"`, `"S4"`, `"S2"`,
#'   or `"all"` for every set that is available.
#' @param quiet If `TRUE`, suppresses progress output and messages.
#' @param .consented Internal. `TRUE` when the caller has already obtained
#'   consent (`.load_bundle()` does so after its own prompt). It skips this
#'   function's prompt and the `allow_download` check. Not for users.
#'
#' @return Invisibly, the paths of the installed files.
#'
#' @examples
#' # Asks first, then downloads about 106 MB per set.
#' # Non-interactively: options(promis2eq5d.allow_download = TRUE)
#' if (interactive()) {
#'   p2e_download("R09_chain", "S3")
#' }
#' @export
p2e_download <- function(method = "R09_chain", set = c("S3", "S1", "S4", "S2", "all"),
                         quiet = FALSE, .consented = FALSE) {
  set <- match.arg(set)
  m <- .p2e_method(method)
  method <- m$id
  wanted <- if (identical(set, "all")) c("S3", "S1", "S4", "S2") else set

  rows <- lapply(wanted, function(s) {
    r <- .manifest_rows(method, s)
    if (nrow(r) == 1L) return(r)
    if (identical(set, "all")) {
      if (!quiet) message("no download is available yet for ", method, " ", s, "; skipping")
      return(NULL)
    }
    stop("no download is available yet for ", method, " ", set, call. = FALSE)
  })
  rows <- Filter(Negate(is.null), rows)
  if (!length(rows)) {
    stop("no download is available yet for ", method, " ", set, call. = FALSE)
  }
  rows <- do.call(rbind, rows)

  if (!isTRUE(.consented)) {
    mb <- sum(as.numeric(rows$bytes), na.rm = TRUE) / 1e6
    if (.is_interactive()) {
      question <- sprintf("Download %d fit file(s) for %s (about %.0f MB) into %s?",
                          nrow(rows), method, mb, file.path(.cache_dir(), "fits"))
      if (!isTRUE(.ask_yes_no(question))) {
        stop("download of ", method, " fits declined", call. = FALSE)
      }
    } else if (!isTRUE(getOption("promis2eq5d.allow_download", FALSE))) {
      stop("download of ", method, " fits needs consent in a non-interactive session; run ",
           "options(promis2eq5d.allow_download = TRUE) and then p2e_download('", method,
           "', set = '", set, "')", call. = FALSE)
    }
  }

  installed <- character()
  for (i in seq_len(nrow(rows))) {
    installed <- c(installed, .download_fit(rows[i, , drop = FALSE], quiet))
  }
  invisible(installed)
}

# The published fit files. Each row: method, set, url, sha256, md5, bytes.
# The R09_chain fits are archived on Zenodo (https://doi.org/10.5281/zenodo.24199380).
.fit_manifest <- function() {
  sets <- c("S1", "S2", "S3", "S4")
  data.frame(
    method = "R09_chain",
    set = sets,
    url = paste0("https://zenodo.org/records/24199380/files/R09_chain_", sets,
                 "_any.rds?download=1"),
    sha256 = c(
      "b8a540f08b3be03ea3fffdbfc6ca0e8b2d710b42ddb4f8533fc0933ab3d31dfa",
      "af03c422b7b640e1ac39100b1ee7e57981107d5cb3c10e62cd464c35877cf696",
      "42d32fbd983278d166f7de8fa87906fc16f3e651631adb9c175a311dbfa2b978",
      "9c9dc7e77bc5ed026acc5108aca644ad072dd84fbf23fd2a26d4cdc5527c805d"),
    md5 = c(
      "2678156da915961ca555c606f9eb84d4",
      "1dcf641e9a0c33ec1df244ffb0be5cdb",
      "c5b579b2b86244bfb725544ee62ee97b",
      "c988e4c312ce548cb415a344d249f7dc"),
    bytes = c(106635044, 106341760, 105946824, 106453656),
    stringsAsFactors = FALSE)
}

# The manifest rows for (method, set): zero or one. Duplicates are an error.
.manifest_rows <- function(method, set) {
  man <- .fit_manifest()
  r <- man[man$method == method & man$set == set, , drop = FALSE]
  if (nrow(r) > 1L) {
    stop("the manifest has duplicate rows for ", method, " ", set, call. = FALSE)
  }
  r
}

.is_interactive <- function() interactive()

.ask_yes_no <- function(question) utils::askYesNo(question, default = FALSE)

# c(sha256 = , md5 = ); sha256 is NA when tools::sha256sum is not available (R < 4.5)
.file_hash <- function(path) {
  sha <- if (exists("sha256sum", envir = asNamespace("tools"), inherits = FALSE)) {
    unname(tools::sha256sum(path))
  } else {
    NA_character_
  }
  c(sha256 = sha, md5 = unname(tools::md5sum(path)))
}

# Compares sha256 when both values are known, else md5
.hash_matches <- function(row, local) {
  has <- function(x) !is.null(x) && !is.na(x) && nzchar(x)
  if (has(row$sha256) && !is.na(local[["sha256"]])) {
    identical(tolower(row$sha256), tolower(local[["sha256"]]))
  } else {
    has(row$md5) && identical(tolower(row$md5), tolower(local[["md5"]]))
  }
}

# One download with utils::download.file(mode = "wb"). Warnings are captured
# and folded into the error when the download fails.
.download_file <- function(url, dest, quiet) {
  warns <- character()
  fail <- function(msg) {
    stop(msg, if (length(warns)) paste0(" [", paste(unique(warns), collapse = "; "), "]"),
         call. = FALSE)
  }
  status <- tryCatch(
    withCallingHandlers(
      utils::download.file(url, dest, mode = "wb", quiet = quiet),
      warning = function(w) {
        warns <<- c(warns, conditionMessage(w))
        invokeRestart("muffleWarning")
      }),
    error = function(e) fail(conditionMessage(e)))
  if (!identical(as.numeric(status), 0)) fail(paste0("exit status ", status))
  invisible(status)
}

# Downloads one manifest row into a temp file, verifies it and installs it in
# the cache. Returns the installed path. A partial or unverified file never
# reaches the final name.
.download_fit <- function(row, quiet) {
  tmp <- tempfile(fileext = ".part")
  on.exit(unlink(tmp), add = TRUE)
  old <- options(timeout = max(getOption("timeout"), 3600))
  on.exit(options(old), add = TRUE)
  tryCatch(
    .download_file(row$url, tmp, quiet),
    error = function(e) {
      stop("download of ", row$method, " ", row$set, " failed: ", conditionMessage(e), call. = FALSE)
    }
  )

  if (!.hash_matches(row, .file_hash(tmp))) {
    stop("checksum mismatch for the downloaded ", row$method, " ", row$set,
         " fit; nothing was installed", call. = FALSE)
  }

  dir <- file.path(.cache_dir(), "fits")
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  final <- file.path(dir, .bundle_file(row$method, row$set, NA))
  staged <- file.path(dir, paste0(".", basename(final), ".", Sys.getpid(), ".part"))
  on.exit(unlink(staged), add = TRUE)

  # Stage in the destination directory (rename, else copy across devices)
  if (!file.rename(tmp, staged) && !file.copy(tmp, staged, overwrite = TRUE)) {
    unlink(staged)
    stop("could not stage the ", row$method, " ", row$set, " fit in ", dir, call. = FALSE)
  }
  if (!file.rename(staged, final)) {
    if (!file.copy(staged, final, overwrite = TRUE)) {
      unlink(final)
      stop("could not install the ", row$method, " ", row$set, " fit at ", final, call. = FALSE)
    }
  }
  if (!quiet) message("installed ", final)
  final
}
