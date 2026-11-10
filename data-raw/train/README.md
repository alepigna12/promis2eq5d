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
`--dirty` become the bundles' `repo_commit` and `repo_dirty`. Folds (`make_cv_folds`, country x util_band, from the row order) are computed once per run
and passed to `.cv_unit()` and `.fit_bundle()` as `fold_id`, so cv units and bundles share them.

`bundle` with no ids rebuilds every bundle (all 23), overwriting existing fits. Keep `--out`
outside the package checkout. Without `--smoke`, `bundle` requires `--commit` and `--dirty 0|1`.
`selected.csv` has a `smoke` column; `bundle` refuses a file made in the other mode. `collect`
shows the target B per bundle and a `short` flag when fewer replicates converged.

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
The scaling also ignores Q (5 in smoke, 50 in the full run, in the out-of-fold simulation),
which adds to R06 and R09 bundle times: give those extra margin.

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

## Running on the JHPCE cluster

Scripts: `data-raw/slurm/submit.sh`, `array.sh`, `cluster.env.example` (adapted from the
iHEA_2027 cluster scripts). No R runs on the login node: `list.sh` (in an `srun` session) lists
unit ids, and `submit.sh` only reads that list and calls `sbatch`.

**Setup (once).**
1. Clone over SSH and switch to the branch:
   `git clone git@github.com:alepigna12/promis2eq5d.git ~/research/promis2eq5d`, then
   `cd ~/research/promis2eq5d && git checkout pkg-train`.
2. Install in an interactive session, not on the login node: `srun --pty bash`, then
   `module load conda_R/4.5`, then `R CMD INSTALL --no-multiarch .` into the user library,
   after installing the package's Suggests there too (`install.packages()` on the list in
   `DESCRIPTION`; some may compile from source).
3. In the same session, run the per-file tests that need no fits (for example
   `Rscript -e 'testthat::test_file("tests/testthat/test-selection.R")'`, and likewise for
   the other files that need no fitted bundles).
4. Copy the licensed T-score table into your own account, outside the repo, and lock it:
   `chmod go-rwx <file>`. Never commit it.
5. `cp data-raw/slurm/cluster.env.example data-raw/slurm/cluster.env` (gitignored) and set
   `TSCORE_TABLE` (and `PROJECT_DIR` if the clone is elsewhere). `OUT` is on fastscratch.
   `TIME_*` are 3 x the projected time of the stage's slowest unit (R06 bundles set
   `TIME_BUNDLE`); `MEM_*` were not measured, so check `sacct` MaxRSS after the cv stage.

   `TIME_BUNDLE` is 48 h because the R06 projection ignores Q; check the shared partition's
   maximum with `sinfo -o "%P %l"` and lower it if needed. `PROJECT_DIR`, `OUT` and
   `TSCORE_TABLE` must not contain commas or spaces (`sbatch --export` is comma-delimited).

**Run order.** No R runs on the login node. `list.sh` (cv and bundle only) runs inside an
`srun --pty` session and writes `$OUT/submit/<stage>.txt`; it loads the installed package, so
the install step above comes first. `submit.sh` runs on the login node (never inside an srun
allocation: `SLURM_*` variables leak through `--export=ALL` and can clash with `--mem`),
reads that list, and calls `sbatch`. Per stage: (a) list inside srun, then exit the session;
(b) the dry-run on the login node; (c) separately, after reading the dry-run, the real
submit. Never chain these with `&&`. A stage starts only after the previous one has finished
(`squeue -u $USER` is empty). cv and bundle:
```
srun --pty --cpus-per-task=1 --mem=4G --time=00:20:00 bash
data-raw/slurm/list.sh cv                  # writes $OUT/submit/cv.txt (153 ids); then exit
DRY_RUN=1 data-raw/slurm/submit.sh cv      # login node
data-raw/slurm/submit.sh cv                # 153 array tasks, one unit each
```
select and collect need no list (a one-element array running `train.R select` / `collect`;
the collect table is in the slurm log):
```
DRY_RUN=1 data-raw/slurm/submit.sh select
data-raw/slurm/submit.sh select
```
Then the same three steps for `bundle` (23 tasks, one bundle each), then the two commands
for `collect`. `SMOKE=1` in front of `list.sh` and `submit.sh` runs the smoke variant
(`list.sh` overwrites `$OUT/submit/<stage>.txt`, so list again when switching). Logs are in
`slurm_logs/` (`p2e_<stage>_<jobid>_<index>.out`); each submission is recorded as
`$OUT/submit/<stamp>_<stage>.txt`. Each task passes `--commit` and `--dirty` from the
cluster checkout (`git rev-parse HEAD`, `git status`), so keep the checkout clean while the
bundle stage runs.

**Afterwards.** fastscratch purges files after 30 days and has no backup: copy `$OUT/fits`
and `$OUT/selected.csv` (and the `slurm_logs`) to your own storage as soon as `collect` is done.
