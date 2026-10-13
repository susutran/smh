#!/bin/bash
#SBATCH --job-name=smh_boltz_prefetch
#SBATCH --account=rgmullen-chbe-eng
#SBATCH --partition=eng-research-gpu
#SBATCH --gres=gpu:1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=02:00:00
#SBATCH --output=logs/boltz-prefetch-%j.out
#SBATCH --error=logs/boltz-prefetch-%j.err
#
# One-off: pull Boltz-2 model weights into a persistent cache.
#
#   sbatch sbatch/5boltz_prefetch.sbatch
#
# Run this ONCE. Every later stage-5 job binds the same cache and downloads
# nothing. Without it, each job re-downloads several GB inside its own
# walltime, which is slow and will eventually fail when a node has no
# outbound network.
#
# Cache location: $SMH_BOLTZ_CACHE, default ~/boltz_cache.
# Home has no purge policy; /scratch is wiped after 30 days of no access, so
# the cache does NOT belong there.

set -euo pipefail

REPO="${SMH_REPO:-${SLURM_SUBMIT_DIR:-$PWD}}"
BOLTZ_SIF="${SMH_BOLTZ_SIF:-$REPO/containers/boltz.sif}"
CACHE="${SMH_BOLTZ_CACHE:-$HOME/boltz_cache}"

[[ -f "$BOLTZ_SIF" ]] || { echo "ERROR: container not found: $BOLTZ_SIF" >&2; exit 1; }

mkdir -p "$CACHE" "$REPO/logs" "$REPO/out/boltz_prefetch"

echo "=== boltz-2 weight prefetch ==="
echo "container: $BOLTZ_SIF"
echo "cache:     $CACHE"
echo "node:      $(hostname)"
echo "started:   $(date -Is)"
echo

echo "--- quota before ---"
du -sh "$CACHE" 2>/dev/null || true
quota -s 2>/dev/null | head -8 || true
echo

# No documented download-only command, so a minimal prediction is what pulls
# the weights. The sequence is a short peptide -- the point is the download,
# not the structure.
WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT
cat > "$WORK/tiny.fasta" << 'FASTA'
>A|protein|empty
MKTAYIAKQRQISFVKSHFSRQ
FASTA

echo "--- running minimal prediction to trigger download ---"
apptainer exec --nv \
    --bind "$CACHE:$CACHE" \
    --bind "$WORK:$WORK" \
    --bind "$REPO:$REPO" \
    "$BOLTZ_SIF" \
    env BOLTZ_CACHE="$CACHE" \
    boltz predict "$WORK/tiny.fasta" \
        --cache "$CACHE" \
        --out_dir "$REPO/out/boltz_prefetch" \
        --use_msa_server \
        --override

echo
echo "--- cache after ---"
du -sh "$CACHE"
find "$CACHE" -maxdepth 2 -type f -size +10M -printf '%10s  %p\n' 2>/dev/null | sort -rn | head -20
echo
echo "finished:  $(date -Is)"
echo
echo "Record the sizes above in the README, and pin the boltz version:"
echo "  apptainer exec $BOLTZ_SIF pip list | grep -Ei 'boltz|torch'"
