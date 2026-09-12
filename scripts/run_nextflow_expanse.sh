#!/usr/bin/env bash
#SBATCH --account=ddp195
#SBATCH --partition=ind-shared
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --job-name=pgs-nextflow-controller

set -euo pipefail

if (( $# == 0 )); then
    echo "usage: run_nextflow_expanse.sh NEXTFLOW_ARGUMENT..." >&2
    exit 2
fi

command -v micromamba >/dev/null
test -x "$HOME/nextflow"

exec micromamba run -n nf_latest \
    env NXF_VER=26.04.6 "$HOME/nextflow" "$@"
