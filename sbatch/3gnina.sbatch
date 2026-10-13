#!/bin/bash
#SBATCH --job-name=smh_3gnina
#SBATCH --account=rgmullen-chbe-eng
#SBATCH --partition=eng-research-gpu
#SBATCH --gres=gpu:1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=logs/3gnina-%A_%a.out
#SBATCH --error=logs/3gnina-%A_%a.err
#
# Stage 3: GPU docking with GNINA.
#
#   sbatch sbatch/3gnina.sbatch <ligands.sdf> <receptor.pdb> <out.sdf>
#
# NOTE: --nv is required for GPU access inside the container. Without it
# GNINA silently falls back to CPU and runs ~50x slower with plausible-
# looking output. The nvidia-smi check below fails loudly instead.

set -euo pipefail

REPO="${SMH_REPO:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
GNINA_SIF="${SMH_GNINA_SIF:-$REPO/containers/gnina.sif}"

if [[ $# -lt 3 ]]; then
    echo "Usage: sbatch $0 <ligands.sdf> <receptor.pdb> <output.sdf> [gnina args]" >&2
    exit 2
fi

LIGANDS="$1"; RECEPTOR="$2"; OUTPUT="$3"; shift 3
for v in LIGANDS RECEPTOR OUTPUT; do
    [[ "${!v}" = /* ]] || declare "$v=$REPO/${!v}"
done

[[ -f "$GNINA_SIF" ]] || { echo "ERROR: $GNINA_SIF not found" >&2; exit 1; }
mkdir -p "$REPO/logs" "$(dirname "$OUTPUT")"

echo "=== smh stage 3: GNINA docking ==="
echo "container: $GNINA_SIF"
echo "ligands:   $LIGANDS"
echo "receptor:  $RECEPTOR"
echo "output:    $OUTPUT"
echo "node:      $(hostname)"
echo

# Fail fast if the GPU is not visible inside the container.
if ! apptainer exec --nv "$GNINA_SIF" nvidia-smi -L; then
    echo "ERROR: no GPU visible inside container. Refusing to run on CPU." >&2
    exit 1
fi
echo

apptainer exec --nv \
    --bind "$REPO:$REPO" \
    "$GNINA_SIF" \
    gnina \
        --receptor "$RECEPTOR" \
        --ligand "$LIGANDS" \
        --out "$OUTPUT" \
        --seed 61453 \
        "$@"

echo
echo "finished:  $(date -Is)"
