This is a new submission.

## Test environments

* Local: Windows 11 x64, R 4.5.2 (2025-10-31 ucrt), x86_64-w64-mingw32.

## R CMD check results

`R CMD check --as-cran` on the built tarball, run with `NOT_CRAN` unset and
`_R_CHECK_CRAN_INCOMING_REMOTE_=false` (no other `_R_CHECK_*` variables set,
so `_R_CHECK_FORCE_SUGGESTS_` has its default): Status OK, 0 errors |
0 warnings | 0 notes. The "New submission" NOTE was not observed locally; it
comes from the incoming checks that CRAN runs, which could not be run here.
The tests gave 1352 passed expectations and 26 skips.

The same check with `_R_CHECK_DEPENDS_ONLY_=true` (recommended packages such
as MASS and Matrix stay available) and the same two settings gives Status OK,
0 errors | 0 warnings | 0 notes. Tests that need a Suggests package are
skipped when it is absent (168 skips in that run).

The installed size is about 3.8 MB (the package ships the fits of three
methods for France, Germany and the UK, and the training data); the source
tarball is 3.3 MB.

## Downloads

The fits of the R09_chain method are large and are not part of the package.
They are archived on Zenodo (doi:10.5281/zenodo.24199380) and are downloaded
by `p2e_download()` only with the user's consent, interactively or by setting
`options(promis2eq5d.allow_download = TRUE)`. Nothing is downloaded during
examples, tests, vignettes or `R CMD check`.
