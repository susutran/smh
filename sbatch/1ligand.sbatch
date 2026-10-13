#!/bin/bash
#SBATCH --job-name=smh_1ligand
#SBATCH --account=rgmullen-chbe-eng
#SBATCH --partition=secondary
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=logs/1ligand-%j.out
#SBATCH --error=logs/1ligand-%j.err
#
# Stage 1: SMILES -> 3D conformers (SDF)
#
# Single run:
#   sbatch sbatch/1ligand.sbatch data/dev_smiles.txt
#
# Array run (10k molecules per task, 20 tasks, 4 concurrent):
#   SMH_CHUNK_SIZE=10000 sbatch --array=0-19%4 sbatch/1ligand.sbatch data/full.txt
#
# Override paths without editing this file:
#   SMH_REPO=/path/to/smh SMH_PREP_SIF=/path/to/prep.sif sbatch ...

set -euo pipefail

REPO="${SMH_REPO:-${SLURM_SUBMIT_DIR:-$PWD}}"
PREP_SIF="${SMH_PREP_SIF:-$REPO/containers/prep.sif}"

if [[ $# -lt 1 ]]; then
    echo "Usage: sbatch $0 <smiles_file> [extra args to 1ligands.py]" >&2
    exit 2
fi

INPUT="$1"; shift
[[ "$INPUT" = /* ]] || INPUT="$REPO/$INPUT"

if [[ ! -f "$PREP_SIF" ]]; then
    echo "ERROR: container not found: $PREP_SIF" >&2
    echo "Build it with: apptainer build --ignore-fakeroot-command $PREP_SIF $REPO/containers/prep.def" >&2
    exit 1
fi

mkdir -p "$REPO/logs"

echo "=== smh stage 1: ligand prep ==="
echo "repo:      $REPO"
echo "container: $PREP_SIF"
echo "input:     $INPUT"
echo "node:      $(hostname)"
echo "started:   $(date -Is)"
echo

apptainer exec \
    --bind "$REPO:$REPO" \
    "$PREP_SIF" \
    python3 "$REPO/src/1ligands.py" "$INPUT" "$@"

echo
echo "finished:  $(date -Is)"
