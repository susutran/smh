#!/usr/bin/env python3
import socket, os, sys
from pathlib import Path

if len(sys.argv) != 2:
    sys.exit("Usage: python 1ligands.py <smiles.txt>")

smiles_path = Path(sys.argv[1])
if not smiles_path.is_file():
    sys.exit(f"Error: file not found: {smiles_path}")
if smiles_path.stat().st_size == 0:
    sys.exit(f"Error: file is empty: {smiles_path}")

print("host:", socket.gethostname())
print("cwd:", os.getcwd())
print("python:", sys.version)
print("SLURM_NTASKS:", os.environ.get("SLURM_NTASKS"))
print("SLURM_JOB_ID:", os.environ.get("SLURM_JOB_ID"))
print("input file:", smiles_path)