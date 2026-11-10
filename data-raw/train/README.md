# Training driver

Fits the shipped bundles of promis2eq5d on `p2e_training` (4,512 rows). It runs on the
installed package (`R CMD INSTALL --no-multiarch .` first) and uses `:::` for internals.
Single process: set `OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1`.

```
Rscript data-raw/train/train.R <stage> --out <dir> --tscore-table <csv> \
    [--smoke] [--commit <sha>] [--dirty 0|1] [unit ids...]
```

The T-score table is the licensed domain/raw/tscore CSV. Pass it only on the command line;
never copy or commit it. Fits (and cv outputs) are built from it and must not be committed
either unless the author decides so.

| stage | what it does |
|---|---|
| `list cv` / `list bundle` | unit ids, one per line (153 cv units, 23 bundle units) |
| `cv [ids]` | one `.cv_unit()` per cell x variant, to `<out>/cv/<id with : as _>.rds`; no ids = every missing unit |
| `select` | fails listing missing cv units, else writes `<out>/selected.csv` (single-variant cells: variant 1) |
| `bundle [ids]` | one `.fit_bundle()` per cell, to `<out>/fits/<method>_<set>_<target>.rds` (xz) and a log in `<out>/bundle_logs/` |
| `collect` | table of every bundle (converged, attempted, failed, variant, xz size, seconds, warnings) and the total size of the CRAN-shipped bundles (all but R09_chain) |

Settings: `TRAIN_SEED = 20261110` (`units.R`); K = 10, Q = 50, B = 200 (R09: 50). Every
write goes to `<file>.tmp` and is renamed. Warnings raised in a unit are collected (not
turned into errors) and stored in the unit's output; the count is printed. `--commit` and
`--dirty` become the bundles' `repo_commit` and `repo_dirty`. `p2e_training` has no ID
column, so folds (`make_cv_folds`, country x util_band) are computed on a copy with an ID
added, and passed to `.fit_bundle()` as `fold_id`; the cv units and bundles share them.

`--smoke`: a stratified 600-row subsample (country x util_band), K = 5, Q = 5, B = 3
(R09: 2), R09 configurations 1-2 only (41 cv units). K is 5 rather than 3 because the D03
lasso variants run `cv.glmnet` on K - 1 inner folds, which needs at least 4; with K = 3 or 4
those units fail (smoke only, full K = 10 is fine). Smoke bundles carry the digest of the
subsample, so they are not valid for `retrain`.

## Time table

Measured on the laptop (smoke, n = 600), scaled to full size (n = 4,512, K = 10).
Scaling: one fit costs about linearly in rows, so x7.5 (4512/600) times the ratio of fits.
cv unit (K fold fits, D03 each with K - 1 inner fits, so K^2 fits): R09 x17 = (10/5) x (4061/480),
D03 x34 (K^2 fits: 100 vs 25, same row factor). Bundle (1 + K + B fits): R06, D03, D12
x176 (211 fits vs 9), R09 x57 (61 vs 8). Superlinear fits (quantreg, mvop) would cost
more; treat these as lower-leaning estimates and set `--time` with margin.

| unit | smoke (s) | full size (est.) |
|---|---|---|
| cv R09 (per config) | 7.3 mean (8 units) | about 2 min (max 3 min) |
| cv D03 (per variant) | 3.0 mean, 8.9 max (33 units) | about 1.7 min mean, 5 min max |
| bundle R06 S1 / S3 / S2 / S4 | 146 / 92 / 100 / 49 | 7.1 h / 4.5 h / 4.9 h / 2.4 h |
| bundle R09 (each) | 5.2 - 5.8 | about 5 min |
| bundle D03 (S1 / S2 / S3 / S4) | 3.5 - 5.2 / 0.9 - 4.8 / 3.3 - 3.6 / 1.4 - 1.6 | 10 - 15 min (S3 10 min, S4 5 min) |
| bundle D12 | 1.3 - 1.8 | about 4 min |

R06 bundles dominate: give them their own long `--time`. Whole smoke run: cv about 3 min,
select seconds, bundle about 8 min.
