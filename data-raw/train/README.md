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
| `list core` / `list merge` | the 16 chunked cells (R06_mvop_dim x 4, D03_clad x 12): `core:<method>:<set>:<target>`, `merge:...` |
| `list boot [--chunk-R06 <n>] [--chunk-D03 <n>]` | the bootstrap chunks of every chunked cell over b = 1..B, `boot:<method>:<set>:<target>:<b1>-<b2>`; chunk sizes default to 20 (R06) and 25 (D03), the last chunk may be shorter (full run: 40 + 96 = 136 chunks) |
| `cv [ids]` | one `.cv_unit()` per cell x variant, to `<out>/cv/<id with : as _>.rds`; no ids = every missing unit |
| `select` | fails listing missing cv units, else writes `<out>/selected.csv` (single-variant cells: variant 1) |
| `bundle [ids]` | one `.fit_bundle()` per cell, to `<out>/fits/<method>_<set>_<target>.rds` (xz) and a log in `<out>/bundle_logs/`; used for R09_chain and D12_klapproth |
| `core [ids]` | chunked cells: `.bundle_core()` (point and fold fits, out-of-fold predictions, the D03 residual pool), to `<out>/core/<bundle file>` (xz) with seconds, warnings, smoke, commit, dirty, variant, B, K, Q |
| `boot [ids]` | `.boot_attempts()` for attempts b1..b2 of a chunked cell (no ids = every chunk of `list boot` with the same chunk options), to `<out>/boot/<bundle file stem>_<b1>-<b2>.rds`; D03 reads the residual pool from its core file (an error when it is missing or from another run) |
| `merge [ids]` | checks that the core and the boot files of the cell exist, have this run's smoke, commit, dirty and variant, and cover 1..B exactly; `.bundle_assemble()` (replacement attempts B + 1, ... as `bundle` makes them) to `<out>/fits/<file>` (xz) and `<out>/bundle_logs/<file>` (seconds = core + chunks + merge, all warnings), so `collect` reads it as a `bundle` output |
| `collect` | table of every bundle (converged, attempted, failed, variant, xz size, seconds, warnings) and the total size of the CRAN-shipped bundles (all but R09_chain) |

Settings: `TRAIN_SEED = 20261110` (`units.R`); K = 10, Q = 50, B = 200 (R09: 50). Every
write goes to `<file>.tmp` and is renamed. Warnings raised in a unit are collected (not
turned into errors) and stored in the unit's output; the count is printed. `--commit` and
`--dirty` become the bundles' `repo_commit` and `repo_dirty`. Folds (`make_cv_folds`, country x util_band, from the row order) are computed once per run
and passed to `.cv_unit()` and `.fit_bundle()` as `fold_id`, so cv units and bundles share them.

`bundle` with no ids rebuilds every bundle (all 23), overwriting existing fits.
The chunked bundles (`CHUNKED_METHODS` in `units.R`: R06_mvop_dim and D03_clad) are fitted as
`core` + `boot` + `merge` instead: the bootstrap attempts are independent (attempt b is fitted
under seed + b on its own resample), so the merged bundle is identical to the one `bundle` makes
in one task (same machine; only `meta$created` differs). `core`, `boot` and `merge` read
`selected.csv`, the training data (with T-scores) and the shared folds like `bundle`, and need
`--commit`/`--dirty` without `--smoke`. Every task recomputes the training data: the core and
boot files hold no training data and no T-score column. Keep `--out`
outside the package checkout. Without `--smoke`, `bundle` requires `--commit` and `--dirty 0|1`.
`selected.csv` has a `smoke` column; `bundle` refuses a file made in the other mode. `collect`
shows the target B per bundle and a `short` flag when fewer replicates converged.

`--smoke`: a stratified 600-row subsample (country x util_band), K = 5, Q = 5, B = 3
(R09: 2), R09 configurations 1-2 only (41 cv units). K is 5 rather than 3 because the D03
lasso variants run `cv.glmnet` on K - 1 inner folds, which needs at least 4; with K = 3 or 4
those units fail (smoke only, full K = 10 is fine). Smoke bundles carry the digest of the
subsample, so they are not valid for `retrain`.

## Time table

Measured on the cluster (full size, cv stage): cv R09 0.5 - 6.5 min per unit (0.85 GB max
RSS); cv D03 16 min to over 3 h per unit (0.5 GB max). These replace the smoke-based cv rows
below, which were far too low for D03 (and R09). The bundle rows below are NOT reliable: they
come from the same scaling and are superseded by the measured `core` times (the core of a
chunked bundle is its point + K fold fits; the core time per fit sets the chunk size).

The original estimates, measured on the laptop (smoke, n = 600), scaled to full size (n = 4,512, K = 10).
Scaling: one fit costs about linearly in rows, so x7.5 (4512/600) times the ratio of fits.
cv unit (K fold fits, D03 each with K - 1 inner fits, so K^2 fits): R09 x17 = (10/5) x (4061/480),
D03 x34 (K^2 fits: 100 vs 25, same row factor). Bundle (1 + K + B fits): R06, D03, D12
x176 (211 fits vs 9), R09 x57 (61 vs 8). Superlinear fits (quantreg, mvop) would cost
more; treat these as lower-leaning estimates and set `--time` with margin.
The scaling also ignores Q (5 in smoke, 50 in the full run, in the out-of-fold simulation),
which adds to R06 and R09 bundle times: give those extra margin.

| unit | smoke (s) | full size (est.) |
|---|---|---|
| cv R09 (per config) | 7.3 mean (8 units) | measured: 0.5 - 6.5 min |
| cv D03 (per variant) | 3.0 mean, 8.9 max (33 units) | measured: 16 min to over 3 h |
| bundle R06 S1 / S3 / S2 / S4 | 146 / 92 / 100 / 49 | unreliable (was 7.1 h / 4.5 h / 4.9 h / 2.4 h); 12 - 23 h or more expected |
| bundle R09 (each) | 5.2 - 5.8 | unreliable (was about 5 min); under 1 h expected |
| bundle D03 (S1 / S2 / S3 / S4) | 3.5 - 5.2 / 0.9 - 4.8 / 3.3 - 3.6 / 1.4 - 1.6 | unreliable (was 10 - 15 min); 1 - 10+ h expected |
| bundle D12 | 1.3 - 1.8 | about 4 min |

R06 and D03 bundles dominate: they are fitted as `core` + `boot` chunks + `merge`. Whole smoke
run: cv about 3 min, select seconds, bundle about 8 min.

## Running on the JHPCE cluster

Scripts: `data-raw/slurm/submit.sh`, `array.sh`, `cluster.env.example` (adapted from the
iHEA_2027 cluster scripts). No R runs on the login node: `list.sh` (in an `srun` session) lists
unit ids, and `submit.sh` only reads that list and calls `sbatch`.

**Setup (once).**
1. Clone over SSH and switch to the branch:
   `git clone git@github.com:alepigna12/promis2eq5d.git ~/research/promis2eq5d`, then
   `cd ~/research/promis2eq5d && git checkout pkg-train`.
2. Install, never on the login node: from the repo root, `mkdir -p slurm_logs`, then
   `sbatch --job-name=p2e_install --partition=shared --cpus-per-task=1 --mem=4G --time=01:00:00 --output=slurm_logs/p2e_install_%j.out data-raw/slurm/install.sh`
   (1 h for the first run, which may compile dependencies; 20 min once they are installed),
   or run `data-raw/slurm/install.sh` inside an `srun --pty` session. It loads `conda_R/4.5`, installs the missing or too-old Imports and
   Suggests into the first library on `.libPaths()` (`~/R/4.5` on JHPCE; some may compile from source), installs this checkout
   with `R CMD INSTALL`, and runs the per-file tests that need no fitted bundles
   (training-data, design, selection). The log is `slurm_logs/p2e_install_<jobid>.out`.
3. Check that the log ends with `done: promis2eq5d ... installed in ...`.
4. Copy the licensed T-score table into your own account, outside the repo, and lock it:
   `chmod go-rwx <file>`. Never commit it.
5. `cp data-raw/slurm/cluster.env.example data-raw/slurm/cluster.env` (gitignored) and set
   `TSCORE_TABLE` (and `PROJECT_DIR` if the clone is elsewhere). `OUT` is on fastscratch.
   `TIME_CORE`/`MEM_CORE`, `TIME_BOOT`/`MEM_BOOT`, `TIME_MERGE`/`MEM_MERGE` and the chunk sizes
   `BOOT_CHUNK_R06`, `BOOT_CHUNK_D03` (empty = 20 and 25) are placeholders: set `TIME_BOOT` and
   the chunk sizes from the measured core times.
   `TIME_*` are 3 x the projected time of the stage's slowest unit (R06 bundles set
   `TIME_BUNDLE`); `MEM_*` were not measured, so check `sacct` MaxRSS after the cv stage.

   `TIME_BUNDLE` is 48 h because the R06 projection ignores Q; check the shared partition's
   maximum with `sinfo -o "%P %l"` and lower it if needed. `PROJECT_DIR`, `OUT` and
   `TSCORE_TABLE` must not contain commas or spaces (`sbatch --export` is comma-delimited).

**Run order.** No R runs on the login node. `list.sh` (cv, bundle, core, boot, merge) runs inside an
`srun --pty` session and writes `$OUT/submit/<stage>.txt`; it loads the installed package, so
the install step above comes first. `submit.sh` runs on the login node (never inside an srun
allocation: `SLURM_*` variables leak through `--export=ALL` and can clash with `--mem`),
reads that list, and calls `sbatch` (one array task per id). Per stage: (a) list inside srun, then exit the session;
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
Then the bundles, in this order (each with the same three steps; the list files can be
filtered with `grep` before the dry-run, since `submit.sh` submits whatever is in
`$OUT/submit/<stage>.txt`):
1. `core` (16 tasks: R06 x 4, D03 x 12). Its seconds per cell (`<out>/core/*.rds`, the slurm
   log) set the chunk sizes and `TIME_BOOT`.
2. `boot`. R06 chunks need no core and may be submitted together with `core` (the one
   exception to waiting for the previous stage); D03 chunks read
   the core's residual pool, so submit them only after the D03 cores have finished. For
   example, inside srun:
   ```
   data-raw/slurm/list.sh boot
   grep ':R06_mvop_dim:' $OUT/submit/boot.txt > $OUT/submit/boot.tmp
   mv $OUT/submit/boot.tmp $OUT/submit/boot.txt      # 40 R06 chunks at the default size
   ```
   and later the same with `':D03_clad:'` (96 chunks), after listing again.
3. `merge` (16 tasks, minutes each) once every core and boot task of the cells has finished.
4. `bundle` only for R09 and D12 (7 tasks): after `list.sh bundle`, keep only those lines:
   ```
   grep -E ':(R09_chain|D12_klapproth):' $OUT/submit/bundle.txt > $OUT/submit/bundle.tmp
   mv $OUT/submit/bundle.tmp $OUT/submit/bundle.txt
   ```

Then the two commands for `collect`. Pull once before the first bundle task, so every core,
boot, merge and bundle task records the same commit (`merge` refuses pieces from another
commit). `SMOKE=1` in front of `list.sh` and `submit.sh` runs the smoke variant
(`list.sh` overwrites `$OUT/submit/<stage>.txt`, so list again when switching). Logs are in
`slurm_logs/` (`p2e_<stage>_<jobid>_<index>.out`); each submission is recorded as
`$OUT/submit/<stamp>_<stage>.txt`. Each task passes `--commit` and `--dirty` from the
cluster checkout (`git rev-parse HEAD`, `git status`), so keep the checkout clean while the
bundle stage runs.

**Afterwards.** fastscratch purges files after 30 days and has no backup: copy `$OUT/fits`
and `$OUT/selected.csv` (and the `slurm_logs`) to your own storage as soon as `collect` is done.
