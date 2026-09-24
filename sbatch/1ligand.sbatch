#!/bin/bash
#SBATCH --time=00:10:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --job-name=smh_1ligand
#SBATCH --account=rgmullen-chbe-eng
#SBATCH --partition=eng-research-gpu
#SBATCH --output=/u/susunt2/smh/logs/1ligand-%j.out
#SBATCH --error=/u/susunt2/smh/logs/1ligand-%j.err

REPO="${SMH_REPO:-/u/susunt2/smh}"
PREP_SIF="${SMH_PREP_SIF:-/scratch/susunt2/prep.sif}"

apptainer exec "$PREP_SIF" python3 "$REPO/src/1ligands.py" "$1"
