#!/bin/bash
# Writes the unit ids of a stage to $OUT/submit/<stage>.txt, one per line:
#   data-raw/slurm/list.sh <cv|bundle|core|boot|merge>
# Run it inside an interactive session (srun --pty --cpus-per-task=1 --mem=4G
# --time=00:20:00 bash), never on the login node: it starts R and loads the installed
# package. Then exit the session and run data-raw/slurm/submit.sh <stage> on the login
# node. SMOKE=1 lists the smoke run's units. Same cluster.env as submit.sh.
# boot: the chunk sizes are BOOT_CHUNK_R06 and BOOT_CHUNK_D03 from cluster.env
# (bootstrap attempts per chunk; empty = the driver's defaults, 20 and 25).
# Adapted from the iHEA_2027 scripts slurm/submit.sh and slurm/array.sh at tag
# freeze-ihea2027-abstract-v3.

set -euo pipefail

stage="${1:-}"
case "$stage" in
  cv|bundle|core|boot|merge) ;;
  *) echo "usage: data-raw/slurm/list.sh <cv|bundle|core|boot|merge>" >&2; exit 1 ;;
esac

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
: "${PROJECT_DIR:?}" "${OUT:?}"

if [ -n "${R_MODULE:-}" ]; then
  # the module's activation scripts read unset variables (JHPCE conda_R:
  # SYS_SYSROOT), so strict mode is off while it loads
  set +u
  module load "$R_MODULE"
  set -u
fi
command -v Rscript >/dev/null 2>&1 || { echo "Rscript not found: set R_MODULE in $env_file" >&2; exit 1; }
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1

list_opt=()
if [ "${SMOKE:-0}" = "1" ]; then list_opt=(--smoke); fi
if [ "$stage" = "boot" ]; then
  if [ -n "${BOOT_CHUNK_R06:-}" ]; then list_opt+=(--chunk-R06 "$BOOT_CHUNK_R06"); fi
  if [ -n "${BOOT_CHUNK_D03:-}" ]; then list_opt+=(--chunk-D03 "$BOOT_CHUNK_D03"); fi
fi

mkdir -p "$OUT/submit"
list="$OUT/submit/${stage}.txt"
cd "$PROJECT_DIR"
Rscript data-raw/train/train.R list "$stage" ${list_opt[@]+"${list_opt[@]}"} \
  | tr -d '\r' | { grep -v '^[[:space:]]*$' || true; } > "$list.tmp"
mv "$list.tmp" "$list"
echo "$list: $(wc -l < "$list" | tr -d ' ') unit(s)"
