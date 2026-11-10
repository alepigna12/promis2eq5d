#!/usr/bin/env bash
#SBATCH --job-name=p2e_install
#SBATCH --partition=shared
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=01:00:00
#SBATCH --output=slurm_logs/p2e_install_%j.out
#
# Installs this checkout of promis2eq5d, and the packages the training run and its
# tests need, into the first library R searches (~/R/4.5 on JHPCE), then runs the per-file tests that
# need no fitted bundles.
#
# Never on the login node. From the repo root, either as a batch job:
#   mkdir -p slurm_logs
#   sbatch --job-name=p2e_install --partition=shared --cpus-per-task=1 --mem=4G #     --time=00:20:00 --output=slurm_logs/p2e_install_%j.out data-raw/slurm/install.sh
# (20 min once the dependencies are installed; the first run, which may compile
# some from source, needs about 1 h) or inside an interactive session:
#   srun --pty --cpus-per-task=1 --mem=4G --time=00:20:00 bash
#   data-raw/slurm/install.sh
#
# R_MODULE defaults to conda_R/4.5 (as in cluster.env.example). Packages go into
# the first library on .libPaths() (on JHPCE ~/R/4.5), where later jobs find them.

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

# Install into the first library R itself searches, so every later R process (the
# R CMD INSTALL test load, the training jobs) finds the packages. On JHPCE the
# conda_R setup puts ~/R/<version> first and leaves R_LIBS_USER off the search
# path, so R_LIBS_USER is not used.
echo "R: $(Rscript -e 'cat(R.version.string)')"
echo "library paths:"
Rscript -e 'cat(paste0("  ", .libPaths()), sep = "\n")'
LIB=$(Rscript -e 'cat(.libPaths()[1])')
RHOME=$(Rscript -e 'cat(R.home())')
case "$LIB" in
  "$RHOME"*)
    echo "error: R's first library ($LIB) is inside the R installation, not a personal" >&2
    echo "       library; create one (e.g. mkdir -p ~/R/4.5) and check .libPaths()" >&2
    exit 1 ;;
esac
[ -w "$LIB" ] || { echo "error: library $LIB is not writable" >&2; exit 1; }
echo "installing into: $LIB"

# Dependencies: Imports and Suggests from DESCRIPTION. Only missing or too-old
# packages are installed; the module may already provide some.
Rscript - "$LIB" <<'EOF'
lib <- commandArgs(trailingOnly = TRUE)[1]
.libPaths(c(lib, .libPaths()))
floors <- c(eq5d = "0.17.0", withr = "0", glmnet = "0", MASS = "0", Matrix = "0",
            ordinal = "0", pbivnorm = "0", quantreg = "0", testthat = "3.2.0",
            tibble = "0", xgboost = "3.2.0")
# Read the installed version from the libraries without loading the package: a
# loaded namespace would keep reporting the old version after an upgrade.
ver <- function(p) {
  v <- suppressWarnings(utils::packageDescription(p, lib.loc = .libPaths(), fields = "Version"))
  if (is.na(v)) NA_character_ else v
}
ok <- function(p) !is.na(ver(p)) && package_version(ver(p)) >= floors[[p]]
need <- names(floors)[!vapply(names(floors), ok, logical(1))]
if (length(need)) {
  message("installing: ", paste(need, collapse = ", "))
  install.packages(need, lib = lib, repos = "https://cloud.r-project.org")
}
still <- names(floors)[!vapply(names(floors), ok, logical(1))]
if (length(still)) stop("still missing or too old: ", paste(still, collapse = ", "))
for (p in names(floors)) message(sprintf("  %-9s %s", p, ver(p)))
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
