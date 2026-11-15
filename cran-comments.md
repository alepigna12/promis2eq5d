This is a new submission.

## Test environments

* Local: Windows 11 x64, R 4.5.2 (2025-10-31 ucrt), x86_64-w64-mingw32.

## R CMD check results

`R CMD check --as-cran` on the built tarball: 0 errors | 0 warnings | 0 notes
locally (with the remote CRAN incoming check disabled; the "New submission"
NOTE appears only in the remote incoming check, which I could not run locally). The same check with
`_R_CHECK_DEPENDS_ONLY_=true` and `_R_CHECK_NO_RECOMMENDED_=true` gives
0 errors | 0 warnings | 0 notes; tests that need a Suggests package are
skipped when it is absent.

The installed size is about 3.8 MB (the package ships the fits of three
methods for France, Germany and the UK, and the training data); the source
tarball is 3.3 MB.

## Downloads

The fits of the R09_chain method are large and are not part of the package.
They are archived on Zenodo (doi:10.5281/zenodo.24199380) and are downloaded
by `p2e_download()` only with the user's consent, interactively or by setting
`options(promis2eq5d.allow_download = TRUE)`. Nothing is downloaded during
examples, tests, vignettes or `R CMD check`.
