# promis2eq5d 0.0.0.9000

* First development version. `predict_eq5d()` maps PROMIS-29 profiles to EQ-5D-5L utilities in mean, individual and cohort modes, with CV+ conformal intervals, for any country in 'eq5d' (response methods R06 and R09; direct methods D03 and D12 with crosswalk or retraining for other countries).
* Helpers: `p2e_columns()`, `p2e_value_set()`, `p2e_tscores()`, `p2e_download()`, and `summary()` and `as.matrix()` methods for results.
* The package ships no PROMIS T-score conversion table.
