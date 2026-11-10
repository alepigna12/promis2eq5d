#!/bin/bash
# Submits one training stage as a SLURM array:
#   data-raw/slurm/submit.sh <cv|select|bundle|collect>
# Adapted from the iHEA_2027 script slurm/submit.sh at tag freeze-ihea2027-abstract-v3.
# Reads data-raw/slurm/cluster.env (PARTITION, MAX_PARALLEL, PROJECT_DIR, OUT,
# TSCORE_TABLE, TIME_<STAGE>, MEM_<STAGE>; optional ACCOUNT, R_MODULE).
# cv and bundle: the unit ids come from `train.R list <stage>` (cheap, no fitting;
# the only R run on the login node) into $OUT/submit/<stamp>_<stage>.txt, one array
# task per id. select and collect: a one-element array whose single task runs the
# stage (no ids); the list file holds the one line "<stage>".
# SMOKE=1 submits the smoke run (--smoke). DRY_RUN=1 prints the sbatch command only
# (the list file is still written, it is needed for the array size).

set -euo pipefail

stage="${1:-}"
case "$stage" in
  cv|select|bundle|collect) ;;
  *) echo "usage: data-raw/slurm/submit.sh <cv|select|bundle|collect>" >&2; exit 1 ;;
esac
STAGE="$(echo "$stage" | tr '[:lower:]' '[:upper:]')"

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
env_file="${CLUSTER_ENV:-$here/cluster.env}"
if [ ! -f "$env_file" ]; then
  echo "missing $env_file: copy data-raw/slurm/cluster.env.example and fill it in" >&2
  exit 1
fi
set -a
# shellcheck disable=SC1090
source <(tr -d '\r' < "$env_file")     # CRLF-safe
set +a
: "${PARTITION:?}" "${MAX_PARALLEL:?}" "${PROJECT_DIR:?}" "${OUT:?}"
time_var="TIME_${STAGE}"; mem_var="MEM_${STAGE}"
time="${!time_var:?$time_var is not set in $env_file}"
mem="${!mem_var:?$mem_var is not set in $env_file}"
if [ "$stage" = "cv" ] || [ "$stage" = "bundle" ]; then
  : "${TSCORE_TABLE:?TSCORE_TABLE is not set in $env_file}"
fi
# --account only when ACCOUNT is set (not needed on JHPCE)
account_opt=()
if [ -n "${ACCOUNT:-}" ]; then account_opt=(--account="$ACCOUNT"); fi
smoke_opt=()
if [ "${SMOKE:-0}" = "1" ]; then smoke_opt=(--smoke); fi

mkdir -p "$PROJECT_DIR/slurm_logs" "$OUT/submit"
stamp="$(date +%Y%m%d-%H%M%S)"
list="$OUT/submit/${stamp}_${stage}.txt"

if [ "$stage" = "cv" ] || [ "$stage" = "bundle" ]; then
  if [ -n "${R_MODULE:-}" ]; then
    # the module's activation scripts read unset variables (JHPCE conda_R:
    # SYS_SYSROOT), so strict mode is off while it loads
    set +u
    module load "$R_MODULE"
    set -u
  fi
  command -v Rscript >/dev/null 2>&1 || { echo "Rscript not found: set R_MODULE in $env_file" >&2; exit 1; }
  # one cheap Rscript call (lists ids, fits nothing); needs the installed package
  (cd "$PROJECT_DIR" && Rscript data-raw/train/train.R list "$stage" ${smoke_opt[@]+"${smoke_opt[@]}"}) \
    | tr -d '\r' | grep -v '^[[:space:]]*$' > "$list.tmp"
  mv "$list.tmp" "$list"
else
  echo "$stage" > "$list"
fi

njobs="$(wc -l < "$list" | tr -d ' ')"
if [ "$njobs" -eq 0 ]; then
  echo "$list: no units, nothing submitted" >&2
  exit 1
fi
if [ "$stage" = "cv" ] || [ "$stage" = "bundle" ]; then
  array="1-${njobs}%${MAX_PARALLEL}"
else
  array="1-1"
fi

cmd=(sbatch --array="$array" --time="$time" --mem="$mem" --cpus-per-task=1
     --partition="$PARTITION" ${account_opt[@]+"${account_opt[@]}"}
     --job-name="p2e_${stage}" --chdir="$PROJECT_DIR"
     --export="ALL,STAGE=${stage},SUBMIT_LIST=${list},PROJECT_DIR=${PROJECT_DIR},OUT=${OUT},TSCORE_TABLE=${TSCORE_TABLE:-},SMOKE=${SMOKE:-0}"
     "$PROJECT_DIR/data-raw/slurm/array.sh")
echo "$list: $njobs unit(s), array $array, time $time, mem $mem"
if [ "${DRY_RUN:-0}" = "1" ]; then
  printf '%q ' "${cmd[@]}"; echo
else
  "${cmd[@]}"
fi
