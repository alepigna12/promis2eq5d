#!/bin/bash
# Submits one training stage as a SLURM array:
#   data-raw/slurm/submit.sh <cv|select|bundle|core|boot|merge|collect>
# Adapted from the iHEA_2027 script slurm/submit.sh at tag freeze-ihea2027-abstract-v3.
# Reads data-raw/slurm/cluster.env (PARTITION, MAX_PARALLEL, PROJECT_DIR, OUT,
# TSCORE_TABLE, TIME_<STAGE>, MEM_<STAGE>; optional ACCOUNT, R_MODULE).
# Run this on the login node, never inside an srun allocation (SLURM_* variables
# would leak through --export=ALL and can clash with --mem). No R runs here.
# cv, bundle, core, boot and merge: the unit ids are read from $OUT/submit/<stage>.txt, made beforehand
# by data-raw/slurm/list.sh in an srun session, and copied to the record
# $OUT/submit/<stamp>_<stage>.txt, one array task per id. select and collect: a
# one-element array whose single task runs the stage (no ids); the record holds
# the one line "<stage>".
# SMOKE=1 submits the smoke run (--smoke); run list.sh with SMOKE=1 too.
# DRY_RUN=1 prints the sbatch command only (the dated record is still written).
# PROJECT_DIR, OUT and TSCORE_TABLE must not contain commas or spaces (--export
# is comma-delimited).

set -euo pipefail

stage="${1:-}"
case "$stage" in
  cv|select|bundle|core|boot|merge|collect) ;;
  *) echo "usage: data-raw/slurm/submit.sh <cv|select|bundle|core|boot|merge|collect>" >&2; exit 1 ;;
esac
# the stages with one array task per unit id (the others run once, without ids)
per_unit=0
case "$stage" in cv|bundle|core|boot|merge) per_unit=1 ;; esac
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
if [ "$per_unit" = "1" ]; then
  : "${TSCORE_TABLE:?TSCORE_TABLE is not set in $env_file}"
fi
# --account only when ACCOUNT is set (not needed on JHPCE)
account_opt=()
if [ -n "${ACCOUNT:-}" ]; then account_opt=(--account="$ACCOUNT"); fi

mkdir -p "$PROJECT_DIR/slurm_logs" "$OUT/submit"
stamp="$(date +%Y%m%d-%H%M%S)"
list="$OUT/submit/${stamp}_${stage}.txt"

if [ "$per_unit" = "1" ]; then
  # the ids list is made beforehand by data-raw/slurm/list.sh inside an srun session
  # (no R runs here); copy it to the dated record of this submission
  src="$OUT/submit/${stage}.txt"
  if [ ! -f "$src" ]; then
    echo "missing $src: make it inside an interactive session, then exit it and rerun:" >&2
    echo "  srun --pty --cpus-per-task=1 --mem=4G --time=00:20:00 bash" >&2
    echo "  data-raw/slurm/list.sh $stage" >&2
    exit 1
  fi
  # grep exits 1 on an empty list, which pipefail would turn into a silent exit
  tr -d '\r' < "$src" | { grep -v '^[[:space:]]*$' || true; } > "$list"
else
  echo "$stage" > "$list"
fi

njobs="$(wc -l < "$list" | tr -d ' ')"
if [ "$njobs" -eq 0 ]; then
  echo "$list: no units, nothing submitted" >&2
  exit 1
fi
if [ "$per_unit" = "1" ]; then
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
