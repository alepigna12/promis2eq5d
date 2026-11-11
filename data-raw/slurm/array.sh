#!/bin/bash
# One SLURM array job of the training driver: array index i runs the unit id on
# line i of the submit list (one unit per task: cv, bundle, core, boot, merge),
# or, for select and collect, the stage itself (no ids). Adapted from the
# iHEA_2027 script slurm/array.sh at tag freeze-ihea2027-abstract-v3. Submitted by data-raw/slurm/submit.sh, which passes
# STAGE, SUBMIT_LIST, PROJECT_DIR, OUT, TSCORE_TABLE and SMOKE (and the variables
# of cluster.env) through --export=ALL.
#SBATCH --cpus-per-task=1
#SBATCH --output=slurm_logs/%x_%A_%a.out

set -euo pipefail

: "${STAGE:?STAGE is not set}"
: "${SUBMIT_LIST:?SUBMIT_LIST is not set}"
: "${PROJECT_DIR:?PROJECT_DIR is not set}"
: "${OUT:?OUT is not set}"
: "${SLURM_ARRAY_TASK_ID:?not running as a SLURM array job}"
# the stages with one unit id per task (they also read the T-score table)
per_unit=0
case "$STAGE" in cv|bundle|core|boot|merge) per_unit=1 ;; esac

cd "$PROJECT_DIR"
if [ -n "${R_MODULE:-}" ]; then
  # the module's activation scripts read unset variables (JHPCE conda_R:
  # SYS_SYSROOT), so strict mode is off while it loads
  set +u
  module load "$R_MODULE"
  set -u
fi
command -v Rscript >/dev/null 2>&1 || { echo "Rscript not found: set R_MODULE in data-raw/slurm/cluster.env" >&2; exit 1; }

# one thread per task, and BLAS must not add threads
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
# temporary files on the user's scratch space when there is one (JHPCE: MYSCRATCH)
export TMPDIR="${MYSCRATCH:-${TMPDIR:-/tmp}}"

# provenance of the bundles: the cluster checkout
commit="$(git -C "$PROJECT_DIR" rev-parse HEAD)"
dirty=0
if [ -n "$(git -C "$PROJECT_DIR" status --porcelain)" ]; then dirty=1; fi

args=("$STAGE" --out "$OUT")
if [ "$per_unit" = "1" ]; then
  : "${TSCORE_TABLE:?TSCORE_TABLE is not set}"
  args+=(--tscore-table "$TSCORE_TABLE")
fi
if [ "${SMOKE:-0}" = "1" ]; then args+=(--smoke); fi
args+=(--commit "$commit" --dirty "$dirty")

ids=()
if [ "$per_unit" = "1" ]; then
  id="$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$SUBMIT_LIST" | tr -d '\r')"
  if [ -z "${id//[[:space:]]/}" ]; then
    echo "array index $SLURM_ARRAY_TASK_ID: no unit on line $SLURM_ARRAY_TASK_ID of $SUBMIT_LIST"
    exit 0
  fi
  ids=("$id")
fi

echo "array index $SLURM_ARRAY_TASK_ID: stage $STAGE ${ids[*]:-} (commit $commit, dirty $dirty)"
Rscript data-raw/train/train.R "${args[@]}" ${ids[@]+"${ids[@]}"}
