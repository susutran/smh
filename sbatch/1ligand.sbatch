#!/bin/bash
#
#SBATCH --time=00:10:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --job-name=susunt2_smh_dev
#SBATCH --account=rgmullen-chbe-eng
#SBATCH --partition=secondary
#SBATCH --output=/u/susunt2/scratch/smh/%j/myjob.o%j
#SBATCH --error=/u/susunt2/scratch/smh/%j/myjob.e%j
#

mkdir -p /u/susunt2/scratch/smh/${SLURM_JOB_ID}
python3 src/1ligands.py "$1"