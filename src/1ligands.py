#!/usr/bin/env python3
"""
Stage 1: ligand preparation.

Reads a whitespace-delimited SMILES file (header row + columns including a
SMILES column and an id column), generates one 3D conformer per molecule with
RDKit ETKDGv3, optimizes with MMFF94s (UFF fallback), and writes an SDF.

Designed for SLURM array jobs: set SLURM_ARRAY_TASK_ID and --chunk-size (or
SMH_CHUNK_SIZE) and each task processes its own slice of the input.

Failures are logged and skipped, never fatal -- a single bad SMILES in a
million-row tranche must not kill the job.
"""

import argparse
import csv
import os
import sys
import time
from pathlib import Path

from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem

# RDKit is noisy on malformed input; we do our own reporting.
RDLogger.DisableLog("rdApp.*")


def parse_args():
    p = argparse.ArgumentParser(description="Prepare 3D ligands from SMILES.")
    p.add_argument("input", type=Path, help="Whitespace-delimited SMILES file")
    p.add_argument("-o", "--output", type=Path, default=None,
                   help="Output SDF (default: <input stem>_prepped[_<task>].sdf)")
    p.add_argument("--smiles-col", default="SMILES",
                   help="Header name of the SMILES column (default: SMILES)")
    p.add_argument("--id-col", default="zincid",
                   help="Header name of the id column (default: zincid)")
    p.add_argument("--chunk-size", type=int,
                   default=int(os.environ.get("SMH_CHUNK_SIZE", 0)),
                   help="Rows per array task. 0 = process whole file.")
    p.add_argument("--task-id", type=int,
                   default=int(os.environ.get("SLURM_ARRAY_TASK_ID", -1)),
                   help="Array task index (0-based). -1 = not an array job.")
    p.add_argument("--max-attempts", type=int, default=3,
                   help="Embedding retries with different seeds (default: 3)")
    p.add_argument("--seed", type=int, default=0xF00D,
                   help="Base random seed -- fixed for reproducibility")
    p.add_argument("--no-optimize", action="store_true",
                   help="Skip forcefield optimization (faster, lower quality)")
    return p.parse_args()


def read_rows(path, smiles_col, id_col):
    """Yield (id, smiles) pairs. Tolerates tab- or space-delimited files."""
    with open(path, newline="") as fh:
        sample = fh.readline()
        fh.seek(0)
        delim = "\t" if "\t" in sample else None
        reader = csv.DictReader(fh, delimiter=delim) if delim else None

        if reader is None:
            # Whitespace-delimited: parse header manually.
            header = fh.readline().split()
            try:
                si = header.index(smiles_col)
                ii = header.index(id_col)
            except ValueError:
                sys.exit(f"Error: columns {smiles_col!r}/{id_col!r} not in header: {header}")
            for lineno, line in enumerate(fh, start=2):
                parts = line.split()
                if len(parts) <= max(si, ii):
                    print(f"  skip line {lineno}: too few columns", file=sys.stderr)
                    continue
                yield parts[ii], parts[si]
        else:
            if smiles_col not in reader.fieldnames or id_col not in reader.fieldnames:
                sys.exit(f"Error: columns {smiles_col!r}/{id_col!r} not in header: "
                         f"{reader.fieldnames}")
            for row in reader:
                yield row[id_col], row[smiles_col]


def embed(mol, max_attempts, base_seed, optimize):
    """Generate and optimize one conformer. Returns None on failure."""
    mol = Chem.AddHs(mol)

    params = AllChem.ETKDGv3()
    params.useSmallRingTorsions = True
    params.enforceChirality = True

    for attempt in range(max_attempts):
        params.randomSeed = base_seed + attempt
        if AllChem.EmbedMolecule(mol, params) == 0:
            break
    else:
        # Last resort: random coords, which usually works where ETKDG fails.
        params.useRandomCoords = True
        params.randomSeed = base_seed
        if AllChem.EmbedMolecule(mol, params) != 0:
            return None

    if optimize:
        try:
            if AllChem.MMFFHasAllMoleculeParams(mol):
                AllChem.MMFFOptimizeMolecule(mol, mmffVariant="MMFF94s", maxIters=500)
            else:
                AllChem.UFFOptimizeMolecule(mol, maxIters=500)
        except Exception:
            # An unoptimized conformer still docks; a crash here does not.
            pass

    return mol


def main():
    args = parse_args()

    if not args.input.is_file():
        sys.exit(f"Error: file not found: {args.input}")
    if args.input.stat().st_size == 0:
        sys.exit(f"Error: file is empty: {args.input}")

    rows = list(read_rows(args.input, args.smiles_col, args.id_col))
    total = len(rows)

    # Slice for array jobs.
    suffix = ""
    if args.chunk_size > 0 and args.task_id >= 0:
        start = args.task_id * args.chunk_size
        end = start + args.chunk_size
        rows = rows[start:end]
        suffix = f"_{args.task_id:05d}"
        if not rows:
            print(f"Task {args.task_id}: no rows in range [{start}:{end}) of {total}. "
                  f"Nothing to do.")
            return 0

    out_path = args.output or args.input.with_name(
        f"{args.input.stem}_prepped{suffix}.sdf")
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fail_path = out_path.with_suffix(".failed.txt")

    print(f"host:      {os.uname().nodename}")
    print(f"job:       {os.environ.get('SLURM_JOB_ID', 'interactive')}"
          f"{f' task {args.task_id}' if args.task_id >= 0 else ''}")
    print(f"input:     {args.input} ({total} rows)")
    print(f"this task: {len(rows)} molecules")
    print(f"output:    {out_path}")
    print(f"rdkit:     {Chem.rdBase.rdkitVersion}", flush=True)

    n_ok = n_parse_fail = n_embed_fail = 0
    t0 = time.time()

    with Chem.SDWriter(str(out_path)) as writer, open(fail_path, "w") as failed:
        failed.write("id\tsmiles\treason\n")

        for i, (mol_id, smiles) in enumerate(rows, start=1):
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                n_parse_fail += 1
                failed.write(f"{mol_id}\t{smiles}\tparse\n")
                continue

            mol3d = embed(mol, args.max_attempts, args.seed, not args.no_optimize)
            if mol3d is None:
                n_embed_fail += 1
                failed.write(f"{mol_id}\t{smiles}\tembed\n")
                continue

            mol3d.SetProp("_Name", mol_id)
            mol3d.SetProp("SMILES", smiles)
            writer.write(mol3d)
            n_ok += 1

            if i % 1000 == 0:
                rate = i / (time.time() - t0)
                print(f"  {i}/{len(rows)}  ({rate:.1f} mol/s)", flush=True)

    elapsed = time.time() - t0
    print(f"\ndone in {elapsed:.1f}s ({n_ok / elapsed:.1f} mol/s)" if elapsed > 0 else "")
    print(f"  succeeded:    {n_ok}")
    print(f"  parse failed: {n_parse_fail}")
    print(f"  embed failed: {n_embed_fail}")

    if n_parse_fail or n_embed_fail:
        print(f"  failures logged to {fail_path}")
    else:
        fail_path.unlink(missing_ok=True)

    # Nonzero exit only if nothing at all worked -- partial failure is normal
    # at scale and should not fail the array task.
    return 0 if n_ok else 1


if __name__ == "__main__":
    sys.exit(main())