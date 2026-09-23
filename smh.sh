#!/bin/bash
if [ ! -f "$1" ]; then
    echo "Error: filepath not found: $1" >&2
    exit 1
fi

sbatch sbatch/1ligand.sbatch "$1"