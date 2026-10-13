#!/bin/bash
#SBATCH --job-name=smh_4pb
#SBATCH --account=rgmullen-chbe-eng
#SBATCH --partition=secondary
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=logs/4pb-%j.out
#SBATCH --error=logs/4pb-%j.err
#
# Stage 4: PoseBusters physical-validity filtering.
#
#   sbatch sbatch/4posebusters.sbatch out/docked.sdf \
#       -p data/ESR1/1XP1_receptor.pdb \
#       -r data/ESR1/1XP1_receptor_ligand.sdf
#
# Without -p, only intramolecular checks run -- receptor clashes go undetected.
# With -r, redock mode adds an RMSD check against the crystal pose, which is
# how you validate the box before trusting a production run.

set -euo pipefail

REPO="${SMH_REPO:-${SLURM_SUBMIT_DIR:-$PWD}}"
PREP_SIF="${SMH_PREP_SIF:-$REPO/containers/prep.sif}"

[[ $# -ge 1 ]] || { echo "Usage: sbatch $0 <poses.sdf> [-p receptor.pdb] [-r ref.sdf]" >&2; exit 2; }

POSES="$1"; shift
[[ "$POSES" = /* ]] || POSES="$REPO/$POSES"
[[ -f "$PREP_SIF" ]] || { echo "ERROR: container not found: $PREP_SIF" >&2; exit 1; }

mkdir -p "$REPO/logs"

echo "=== smh stage 4: PoseBusters ==="
echo "container: $PREP_SIF"
echo "poses:     $POSES"
echo "node:      $(hostname)"
echo

apptainer exec --bind "$REPO:$REPO" "$PREP_SIF" \
    python3 "$REPO/src/4posebusters.py" "$POSES" "$@"

echo
echo "finished:  $(date -Is)"
