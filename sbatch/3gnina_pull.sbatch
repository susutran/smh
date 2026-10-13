#!/bin/bash
#SBATCH --account=rgmullen-chbe-eng
#SBATCH --partition=secondary
#SBATCH --cpus-per-task=8
#SBATCH --mem=256G
#SBATCH --time=03:00:00
#SBATCH --job-name=gnina-pull
#SBATCH --output=/u/susunt2/smh/pull-%j.log

export APPTAINER_CACHEDIR=/scratch/susunt2/.apptainer
export APPTAINER_TMPDIR=/tmp/susunt2-build
mkdir -p $APPTAINER_TMPDIR

apptainer pull /scratch/susunt2/gnina.sif docker://gnina/gnina:v1.3.1
