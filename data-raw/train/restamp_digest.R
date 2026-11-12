# Re-stamps the training-data digest of fit bundles from the cluster run.
#
# Why this exists: the bundles of the 2026-11-12 cluster run carry a
# meta$train_digest computed by the first .data_digest(), the md5 of a gzip
# saveRDS(version = 3) file. That depends on the writing R version (the file
# header) and the platform (gzip OS byte), so it did not match on other
# machines. .data_digest() now hashes the header-free serialization payload.
# This script recomputes the old digest the way the cluster did, checks it
# against each bundle, and replaces only meta$train_digest with the new one.
#
# Usage: Rscript data-raw/train/restamp_digest.R <in_dir> <out_dir>
# Run from the repository root; it loads the package source with pkgload::load_all().

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("usage: Rscript data-raw/train/restamp_digest.R <in_dir> <out_dir>")
in_dir <- args[[1]]; out_dir <- args[[2]]
suppressMessages(pkgload::load_all(".", quiet = TRUE))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

train <- as.data.frame(promis2eq5d:::.training_data())

# The cluster's digest: saveRDS(version = 3), gzip, with the writer R version
# of the cluster's R and the gzip OS byte of Linux (3)
old_digest <- function(d, r_version) {
  m <- regmatches(r_version, regexec("^R version ([0-9]+)\\.([0-9]+)\\.([0-9]+)", r_version))[[1]]
  if (length(m) != 4L) stop("cannot parse r_version: ", r_version)
  w <- as.integer(m[2:4])
  d <- d[promis2eq5d:::.digest_cols()]
  rownames(d) <- NULL
  s <- serialize(d, NULL, version = 3)
  s[8:10] <- as.raw(w)
  # the gzip file of that stream, then the OS byte of the gzip header
  f <- tempfile(fileext = ".rds"); on.exit(unlink(f), add = TRUE)
  con <- gzfile(f, "wb"); writeBin(s, con); close(con)
  r <- readBin(f, "raw", file.size(f))
  r[10] <- as.raw(3L)
  writeBin(r, f)
  unname(tools::md5sum(f))
}

new_digest <- promis2eq5d:::.data_digest(train)
files <- list.files(in_dir, pattern = "\\.rds$", full.names = TRUE)
if (!length(files)) stop("no .rds files in ", in_dir)
for (p in files) {
  b <- readRDS(p)
  od <- old_digest(train, b$meta$r_version)
  if (!identical(od, b$meta$train_digest)) {
    stop(basename(p), ": the old-style digest of the training data (", od,
         ") does not equal the bundle's (", b$meta$train_digest, ")")
  }
  if (!identical(as.integer(b$meta$n), nrow(train))) stop(basename(p), ": meta$n differs from the training data")
  b2 <- b
  b2$meta$train_digest <- new_digest
  out <- file.path(out_dir, basename(p))
  tmp <- paste0(out, ".tmp")
  saveRDS(b2, tmp, compress = "xz")
  if (!file.rename(tmp, out)) stop("could not rename ", tmp, " to ", out)
  chk <- readRDS(out)
  expected <- b; expected$meta$train_digest <- new_digest
  same <- identical(chk, expected)
  cat(sprintf("%s old=%s new=%s identical=%s\n", basename(p), od, new_digest, same))
  if (!same) stop(basename(p), ": re-read output differs from the input beyond meta$train_digest")
}
