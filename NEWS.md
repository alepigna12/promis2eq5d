# promis2eq5d 0.0.0.9000

* First development version. `predict_eq5d()` maps PROMIS-29 profiles to EQ-5D-5L utilities in mean, individual and cohort modes, with CV+ conformal intervals, for any country in 'eq5d' (response methods R06 and R09; direct methods D03 and D12 with crosswalk or retraining for other countries).
* Helpers: `p2e_columns()`, `p2e_value_set()`, `p2e_tscores()`, `p2e_download()`, and `summary()` and `as.matrix()` methods for results.
* The package ships no PROMIS T-score conversion table.
* `p2e_training`: the CC0 training data (4,512 respondents; PROMIS-29 items, age, sex, country and the EQ-5D-5L responses, no T-scores) are shipped with the package and used for retraining; the age warning now uses the age range of these data.
* The fits of the response method R06_mvop_dim and the direct methods D03_clad and D12_klapproth (France, Germany, UK) are shipped with the package (`inst/fits`), so `predict_eq5d()` runs for these methods without a download; the training-data digest stored in each fit is now portable across R versions and platforms.
