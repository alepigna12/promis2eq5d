# promis2eq5d 0.1.0

First release.

* `predict_eq5d()` maps PROMIS-29 profiles to EQ-5D-5L utilities for any country value set in 'eq5d', with the response methods `R06_mvop_dim` and `R09_chain` and the direct methods `D03_clad` and `D12_klapproth`. Other countries are served by crosswalk or by retraining on the shipped training data.
* Three modes: `"mean"` (expected utility, optional CV+ conformal interval and most likely dimension levels), `"individual"` (draws per row) and `"cohort"` (group mean utility per draw, aligned across groups and time points for QALY calculations), with `id` and `within_id` for several rows per person.
* Helpers: `p2e_columns()`, `p2e_value_set()`, `p2e_tscores()`, `p2e_download()`, and `summary()` and `as.matrix()` methods for results. The package ships no PROMIS T-score conversion table.
* The fits of `R06_mvop_dim`, `D03_clad` and `D12_klapproth` (France, Germany, UK) are shipped with the package. The R09_chain fits are archived on Zenodo (doi:10.5281/zenodo.24199380) and installed with `p2e_download()`, which asks for consent (or needs `options(promis2eq5d.allow_download = TRUE)`) and checks the SHA-256 checksum.
* Data: `p2e_training` (the CC0 training data, 4,512 respondents) and `p2e_example_trial` (a synthetic two-arm trial).
* Three vignettes: getting started, methods, and uncertainty for cost-effectiveness.
* `citation("promis2eq5d")` gives the package, the R09 fits and the training data.
