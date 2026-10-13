#!/bin/bash
#SBATCH --job-name=smh_2protein
#SBATCH --account=rgmullen-chbe-eng
#SBATCH --partition=secondary
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:20:00
#SBATCH --output=logs/2protein-%j.out
#SBATCH --error=logs/2protein-%j.err
#
# Stage 2: receptor cleanup + docking box.
#
# Run --list-het FIRST, interactively, and look at the output before
# committing to a box. This stage takes seconds; getting it wrong invalidates
# every downstream result.
#
#   apptainer exec containers/prep.sif python3 src/2protein.py \
#       data/ESR1/1XP1.cif --list-het
#
# Then:
#   sbatch sbatch/2protein.sbatch data/ESR1/1XP1.cif -o data/ESR1/1XP1_receptor

set -euo pipefail

REPO="${SMH_REPO:-${SLURM_SUBMIT_DIR:-$PWD}}"
PREP_SIF="${SMH_PREP_SIF:-$REPO/containers/prep.sif}"

[[ $# -ge 1 ]] || { echo "Usage: sbatch $0 <structure.cif|pdb> [args]" >&2; exit 2; }

INPUT="$1"; shift
[[ "$INPUT" = /* ]] || INPUT="$REPO/$INPUT"
[[ -f "$PREP_SIF" ]] || { echo "ERROR: container not found: $PREP_SIF" >&2; exit 1; }

mkdir -p "$REPO/logs"

echo "=== smh stage 2: receptor prep ==="
echo "container: $PREP_SIF"
echo "input:     $INPUT"
echo "node:      $(hostname)"
echo

apptainer exec --bind "$REPO:$REPO" "$PREP_SIF" \
    python3 "$REPO/src/2protein.py" "$INPUT" "$@"

echo
echo "finished:  $(date -Is)"
