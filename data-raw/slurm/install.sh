#!/usr/bin/env bash
#SBATCH --job-name=p2e_install
#SBATCH --partition=shared
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --output=slurm_logs/p2e_install_%j.out
#
# Installs this checkout of promis2eq5d, and the packages the training run and its
# tests need, into the R user library on JHPCE, then runs the per-file tests that
# need no fitted bundles.
#
# Never on the login node. From the repo root, either as a batch job:
#   mkdir -p slurm_logs
#   sbatch data-raw/slurm/install.sh
# or inside an interactive session:
#   srun --pty --cpus-per-task=1 --mem=8G --time=02:00:00 bash
#   data-raw/slurm/install.sh
#
# R_MODULE defaults to conda_R/4.5 (as in cluster.env.example).

set -euo pipefail

if [ -z "${SLURM_JOB_ID:-}" ]; then
  echo "error: run inside a Slurm job (sbatch data-raw/slurm/install.sh, or inside" >&2
  echo "       srun --pty bash), never on the login node" >&2
  exit 1
fi
if [ ! -f DESCRIPTION ] || ! grep -q '^Package: promis2eq5d$' DESCRIPTION; then
  echo "error: run from the root of the promis2eq5d checkout" >&2
  exit 1
fi

R_MODULE=${R_MODULE:-conda_R/4.5}
set +u
module load "$R_MODULE"
set -u
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1

# The user library, created if missing (install.packages() would otherwise ask)
LIB=$(Rscript -e 'cat(path.expand(Sys.getenv("R_LIBS_USER")))')
mkdir -p "$LIB"
echo "R: $(Rscript -e 'cat(R.version.string)')"
echo "user library: $LIB"

# Dependencies: Imports and Suggests from DESCRIPTION. Only missing or too-old
# packages are installed; the module may already provide some.
Rscript - "$LIB" <<'EOF'
lib <- commandArgs(trailingOnly = TRUE)[1]
.libPaths(c(lib, .libPaths()))
floors <- c(eq5d = "0.17.0", withr = "0", glmnet = "0", MASS = "0", Matrix = "0",
            ordinal = "0", pbivnorm = "0", quantreg = "0", testthat = "3.2.0",
            tibble = "0", xgboost = "3.2.0")
ok <- function(p) requireNamespace(p, quietly = TRUE) &&
  utils::packageVersion(p) >= floors[[p]]
need <- names(floors)[!vapply(names(floors), ok, logical(1))]
if (length(need)) {
  message("installing: ", paste(need, collapse = ", "))
  install.packages(need, lib = lib, repos = "https://cloud.r-project.org")
}
still <- names(floors)[!vapply(names(floors), ok, logical(1))]
if (length(still)) stop("still missing or too old: ", paste(still, collapse = ", "))
for (p in names(floors)) message(sprintf("  %-9s %s", p, utils::packageVersion(p)))
EOF

R CMD INSTALL --no-multiarch --library="$LIB" .

# Per-file tests that need no fitted bundles
Rscript - "$LIB" <<'EOF'
lib <- commandArgs(trailingOnly = TRUE)[1]
.libPaths(c(lib, .libPaths()))
for (f in c("test-training-data.R", "test-design.R", "test-selection.R")) {
  message("== ", f)
  testthat::test_file(file.path("tests", "testthat", f), package = "promis2eq5d",
                      load_package = "installed", stop_on_failure = TRUE)
}
EOF

echo "done: promis2eq5d $(Rscript -e 'cat(format(utils::packageVersion("promis2eq5d")))') installed in $LIB"
