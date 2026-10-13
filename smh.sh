#!/bin/bash
#
# smh pipeline driver.
#
# Submits every stage as a SLURM job with afterok dependencies, so the whole
# screen runs unattended and a failed stage cancels everything downstream
# instead of running on a missing or half-written file.
#
#   ./smh.sh data/dev_smiles.txt data/ESR1/1XP1.cif
#   ./smh.sh data/dev_smiles.txt data/ESR1/1XP1.cif --dry-run
#   ./smh.sh data/dev_smiles.txt data/ESR1/1XP1.cif --boltz
#
# Every run gets its own directory under out/, holding the outputs, the job
# ids, and a manifest recording the git commit, container versions and inputs
# that produced them. A result you cannot trace back to a commit and a
# container is not reproducible, whatever the README says.
#
# Stage graph:
#
#     0 filter ──► 1 ligands ──┐
#                              ├──► 3 gnina ──► 4 posebusters ──► [5 boltz]
#     2 receptor ──────────────┘
#
# Stage 2 runs in parallel with 0/1 -- the receptor does not depend on the
# ligands. It is skipped if the prepared receptor already exists, since it is
# per-target, not per-screen.
#
# KNOWN LIMIT: single job per stage, no SLURM arrays. Fine up to ~10^4
# molecules. Beyond that, stages 0 and 1 need --array plus a merge step
# before stage 3. See README.

set -euo pipefail

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$REPO"

usage() {
    cat >&2 << 'USAGE'
Usage: ./smh.sh <smiles_file> <receptor.cif|pdb> [options]

Options:
  --boltz              also run stage 5 (Boltz-2 co-folding)
  --skip-filter        skip stage 0, treat the input as already filtered
  --filter-config F    JSON config for stage 0 thresholds
  --autobox-ligand L   HET residue name for the docking box (default: largest)
  --num-modes N        poses per ligand from GNINA (default: 9)
  --run-name NAME      output directory name (default: timestamp)
  --dry-run            print the submissions without running them
  -h, --help           this message
USAGE
    exit 2
}

[[ $# -ge 2 ]] || usage

SMILES="$1"; RECEPTOR="$2"; shift 2

BOLTZ=0
SKIP_FILTER=0
FILTER_CONFIG=""
AUTOBOX_LIGAND=""
NUM_MODES=9
RUN_NAME=""
DRY_RUN=0

while [[ $# -gt 0 ]]; do
    case "$1" in
        --boltz)           BOLTZ=1; shift ;;
        --skip-filter)     SKIP_FILTER=1; shift ;;
        --filter-config)   FILTER_CONFIG="$2"; shift 2 ;;
        --autobox-ligand)  AUTOBOX_LIGAND="$2"; shift 2 ;;
        --num-modes)       NUM_MODES="$2"; shift 2 ;;
        --run-name)        RUN_NAME="$2"; shift 2 ;;
        --dry-run)         DRY_RUN=1; shift ;;
        -h|--help)         usage ;;
        *) echo "Unknown option: $1" >&2; usage ;;
    esac
done

[[ "$SMILES"   = /* ]] || SMILES="$REPO/$SMILES"
[[ "$RECEPTOR" = /* ]] || RECEPTOR="$REPO/$RECEPTOR"

for f in "$SMILES" "$RECEPTOR"; do
    [[ -f "$f" ]] || { echo "ERROR: not found: $f" >&2; exit 1; }
done

PREP_SIF="${SMH_PREP_SIF:-$REPO/containers/prep.sif}"
GNINA_SIF="${SMH_GNINA_SIF:-$REPO/containers/gnina.sif}"
BOLTZ_SIF="${SMH_BOLTZ_SIF:-$REPO/containers/boltz.sif}"

for sif in "$PREP_SIF" "$GNINA_SIF"; do
    [[ -f "$sif" ]] || { echo "ERROR: container not found: $sif" >&2; exit 1; }
done
if [[ $BOLTZ -eq 1 && ! -f "$BOLTZ_SIF" ]]; then
    echo "ERROR: --boltz given but $BOLTZ_SIF not found" >&2; exit 1
fi

RUN_NAME="${RUN_NAME:-$(date +%Y%m%d-%H%M%S)}"
RUN="$REPO/out/$RUN_NAME"
mkdir -p "$RUN" "$REPO/logs"

RECEPTOR_STEM="${RECEPTOR%.*}_receptor"
RECEPTOR_PDB="${RECEPTOR_STEM}.pdb"
RECEPTOR_LIG="${RECEPTOR_STEM}_ligand.pdb"

SMILES_STEM="$(basename "${SMILES%.*}")"
FILTERED="$RUN/${SMILES_STEM}_filtered.txt"
LIGANDS="$RUN/${SMILES_STEM}_prepped.sdf"
DOCKED="$RUN/docked.sdf"

# --- provenance ----------------------------------------------------------
MANIFEST="$RUN/manifest.txt"
{
    echo "run:        $RUN_NAME"
    echo "started:    $(date -Is)"
    echo "user:       $USER"
    echo "repo:       $REPO"
    echo "git commit: $(git -C "$REPO" rev-parse HEAD 2>/dev/null || echo 'not a git repo')"
    echo "git status: $(git -C "$REPO" status --porcelain 2>/dev/null | wc -l) uncommitted file(s)"
    echo
    echo "smiles:     $SMILES"
    echo "receptor:   $RECEPTOR"
    echo "autobox:    ${AUTOBOX_LIGAND:-auto (largest non-junk HET)}"
    echo "num_modes:  $NUM_MODES"
    echo "filter:     $([[ $SKIP_FILTER -eq 1 ]] && echo skipped || echo "${FILTER_CONFIG:-defaults}")"
    echo "boltz:      $([[ $BOLTZ -eq 1 ]] && echo yes || echo no)"
    echo
    echo "containers:"
    for sif in "$PREP_SIF" "$GNINA_SIF" $([[ $BOLTZ -eq 1 ]] && echo "$BOLTZ_SIF"); do
        echo "  $sif  $(stat -c '%s bytes  mtime %y' "$sif")"
    done
} > "$MANIFEST"

echo "=== smh pipeline ==="
cat "$MANIFEST"
echo

submit() {
    # submit <description> <dependency-or-empty> <sbatch args...>
    local desc="$1"; shift
    local dep="$1"; shift
    local depflag=()
    [[ -n "$dep" ]] && depflag=(--dependency="afterok:$dep")

    if [[ $DRY_RUN -eq 1 ]]; then
        echo "[dry-run] $desc" >&2
        echo "          sbatch ${depflag[*]} $*" >&2
        echo "000000"
        return
    fi

    local jid
    jid=$(sbatch --parsable "${depflag[@]}" "$@")
    echo "  $desc -> job $jid${dep:+  (after $dep)}" >&2
    echo "$jid"
}

# --- stage 2: receptor (independent of ligands) --------------------------
JOB_RECEPTOR=""
if [[ -f "$RECEPTOR_PDB" && -f "$RECEPTOR_LIG" ]]; then
    echo "  stage 2 receptor -> already prepared, skipping ($RECEPTOR_PDB)" >&2
else
    args=("$REPO/sbatch/2protein.sbatch" "$RECEPTOR" -o "$RECEPTOR_STEM")
    [[ -n "$AUTOBOX_LIGAND" ]] && args+=(--autobox-ligand "$AUTOBOX_LIGAND")
    JOB_RECEPTOR=$(submit "stage 2 receptor" "" "${args[@]}")
fi

# --- stage 0: filters ----------------------------------------------------
JOB_PREV=""
LIGAND_INPUT="$SMILES"
if [[ $SKIP_FILTER -eq 0 ]]; then
    args=("$REPO/sbatch/0filter.sbatch" "$SMILES" -o "$FILTERED")
    [[ -n "$FILTER_CONFIG" ]] && args+=(--config "$FILTER_CONFIG")
    JOB_PREV=$(submit "stage 0 filter" "" "${args[@]}")
    LIGAND_INPUT="$FILTERED"
fi

# --- stage 1: ligand prep ------------------------------------------------
# Stage 0 writes a 2-column file; tell stage 1 what the columns are called.
JOB_LIGANDS=$(submit "stage 1 ligands" "$JOB_PREV" \
    "$REPO/sbatch/1ligand.sbatch" "$LIGAND_INPUT" -o "$LIGANDS" \
    --id-col zincid --smiles-col SMILES)

# --- stage 3: docking ----------------------------------------------------
DEP3="$JOB_LIGANDS"
[[ -n "$JOB_RECEPTOR" ]] && DEP3="$JOB_LIGANDS:$JOB_RECEPTOR"

JOB_DOCK=$(submit "stage 3 gnina" "$DEP3" \
    "$REPO/sbatch/3gnina.sbatch" "$LIGANDS" "$RECEPTOR_PDB" "$DOCKED" \
    --autobox_ligand "$RECEPTOR_LIG" --num_modes "$NUM_MODES")

# --- stage 4: posebusters ------------------------------------------------
JOB_PB=$(submit "stage 4 posebusters" "$JOB_DOCK" \
    "$REPO/sbatch/4posebusters.sbatch" "$DOCKED" \
    -p "$RECEPTOR_PDB" -o "$RUN/pb")

# --- stage 5: boltz (optional) -------------------------------------------
JOB_BOLTZ=""
if [[ $BOLTZ -eq 1 ]]; then
    if [[ -f "$REPO/sbatch/5boltz.sbatch" ]]; then
        JOB_BOLTZ=$(submit "stage 5 boltz" "$JOB_PB" \
            "$REPO/sbatch/5boltz.sbatch" "$RUN/pb_passed.sdf" "$RECEPTOR_PDB" \
            -o "$RUN/boltz")
    else
        echo "  stage 5 boltz -> sbatch/5boltz.sbatch not written yet, skipping" >&2
    fi
fi

{
    echo
    echo "jobs:"
    [[ -n "$JOB_RECEPTOR" ]] && echo "  stage 2 receptor    $JOB_RECEPTOR"
    [[ -n "$JOB_PREV"     ]] && echo "  stage 0 filter      $JOB_PREV"
    echo "  stage 1 ligands     $JOB_LIGANDS"
    echo "  stage 3 gnina       $JOB_DOCK"
    echo "  stage 4 posebusters $JOB_PB"
    [[ -n "$JOB_BOLTZ"    ]] && echo "  stage 5 boltz       $JOB_BOLTZ"
} | tee -a "$MANIFEST"

cat << EOF

Submitted. Outputs -> $RUN

  watch:  squeue -u $USER
  logs:   tail -f $REPO/logs/*-${JOB_LIGANDS}.out
  cancel: scancel ${JOB_PREV:-}${JOB_PREV:+ }$JOB_LIGANDS $JOB_DOCK $JOB_PB ${JOB_RECEPTOR:-} ${JOB_BOLTZ:-}

A failed stage cancels everything downstream (afterok). If a job vanishes
from squeue without output, check sacct for its state -- a dependency that
can never be satisfied shows as DependencyNeverSatisfied, not as an error.
EOF