#!/bin/bash
#SBATCH --job-name=smh_0filter
#SBATCH --account=rgmullen-chbe-eng
#SBATCH --partition=secondary
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=logs/0filter-%j.out
#SBATCH --error=logs/0filter-%j.err
#
# Stage 0: cheap 2D filters. Runs before conformer generation.
#
#   sbatch sbatch/0filter.sbatch data/dev_smiles.txt
#   sbatch sbatch/0filter.sbatch data/full.txt --config configs/leadlike.json
#
# Array over a large tranche:
#   SMH_CHUNK_SIZE=200000 sbatch --array=0-49%8 sbatch/0filter.sbatch data/full.txt

set -euo pipefail

REPO="${SMH_REPO:-${SLURM_SUBMIT_DIR:-$PWD}}"
PREP_SIF="${SMH_PREP_SIF:-$REPO/containers/prep.sif}"

[[ $# -ge 1 ]] || { echo "Usage: sbatch $0 <smiles_file> [args]" >&2; exit 2; }

INPUT="$1"; shift
[[ "$INPUT" = /* ]] || INPUT="$REPO/$INPUT"
[[ -f "$PREP_SIF" ]] || { echo "ERROR: container not found: $PREP_SIF" >&2; exit 1; }

mkdir -p "$REPO/logs"

echo "=== smh stage 0: 2D filters ==="
echo "container: $PREP_SIF"
echo "input:     $INPUT"
echo "node:      $(hostname)"
echo "started:   $(date -Is)"
echo

apptainer exec --bind "$REPO:$REPO" "$PREP_SIF" \
    python3 "$REPO/src/0filter.py" "$INPUT" "$@"

echo
echo "finished:  $(date -Is)"
