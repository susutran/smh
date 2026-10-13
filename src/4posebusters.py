#!/usr/bin/env python3
"""
Stage 4: PoseBusters physical-plausibility filtering.

Docking scores do not check whether a pose is physically possible. GNINA will
happily return a high-CNN-score pose with a distorted ring, an impossible bond
angle, or atoms overlapping the protein. PoseBusters catches those.

Reads the docked SDF, runs PoseBusters in redock mode when a reference ligand
is supplied (adds an RMSD check against the crystal pose) or dock mode
otherwise, and writes:

  <out>_passed.sdf   poses that pass every check
  <out>_report.csv   per-pose, per-check results
  <out>_summary.txt  failure counts by check

Exit code is 0 as long as the run completed, even if nothing passed -- a 0%
pass rate is a real and informative result, usually about the receptor or box
rather than the ligands.
"""

import argparse
import sys
from pathlib import Path

import pandas as pd
from rdkit import Chem, RDLogger

RDLogger.DisableLog("rdApp.*")

try:
    from posebusters import PoseBusters
except ImportError:
    sys.exit("Error: posebusters not installed in this container.\n"
             "Rebuild prep.sif with posebusters in containers/prep.def.")


def parse_args():
    p = argparse.ArgumentParser(description="Filter docked poses by physical validity.")
    p.add_argument("poses", type=Path, help="Docked poses (SDF)")
    p.add_argument("-p", "--protein", type=Path, default=None,
                   help="Receptor PDB. Required for protein-ligand checks; "
                        "without it only intramolecular checks run.")
    p.add_argument("-r", "--reference", type=Path, default=None,
                   help="Crystal ligand (SDF/PDB). Enables redock mode + RMSD.")
    p.add_argument("-o", "--output", type=Path, default=None,
                   help="Output stem (default: <poses stem>_pb)")
    p.add_argument("--rmsd-threshold", type=float, default=2.0,
                   help="RMSD cutoff in redock mode (default: 2.0 A)")
    p.add_argument("--top-n", type=int, default=0,
                   help="Only check the first N poses per molecule (0 = all)")
    p.add_argument("--keep-failed", action="store_true",
                   help="Also write the failing poses to <out>_failed.sdf")
    return p.parse_args()


def main():
    args = parse_args()

    if not args.poses.is_file():
        sys.exit(f"Error: file not found: {args.poses}")
    for f in (args.protein, args.reference):
        if f is not None and not f.is_file():
            sys.exit(f"Error: file not found: {f}")

    # Mode selection drives which checks run. Be explicit about it -- "dock"
    # silently skips the RMSD check, which is the one people assume they got.
    if args.reference:
        mode = "redock"
    elif args.protein:
        mode = "dock"
    else:
        mode = "mol"
        print("WARNING: no --protein given. Running intramolecular checks only; "
              "clashes with the receptor will NOT be detected.", file=sys.stderr)

    out_stem = args.output or args.poses.with_name(args.poses.stem + "_pb")
    out_stem = Path(out_stem)
    out_stem.parent.mkdir(parents=True, exist_ok=True)

    supplier = Chem.SDMolSupplier(str(args.poses), removeHs=False, sanitize=True)
    mols, skipped = [], 0
    per_mol = {}
    for mol in supplier:
        if mol is None:
            skipped += 1
            continue
        if args.top_n:
            name = mol.GetProp("_Name") if mol.HasProp("_Name") else ""
            per_mol[name] = per_mol.get(name, 0) + 1
            if per_mol[name] > args.top_n:
                continue
        mols.append(mol)

    if not mols:
        sys.exit(f"Error: no readable poses in {args.poses} ({skipped} unparseable)")

    print(f"poses:     {args.poses} ({len(mols)} loaded, {skipped} unreadable)")
    print(f"mode:      {mode}")
    print(f"protein:   {args.protein or 'none'}")
    print(f"reference: {args.reference or 'none'}", flush=True)

    buster = PoseBusters(config=mode)
    df = buster.bust(mols,
                     mol_true=str(args.reference) if args.reference else None,
                     mol_cond=str(args.protein) if args.protein else None,
                     full_report=True)

    # Boolean check columns only; the rest are metadata/measurements.
    check_cols = [c for c in df.columns if df[c].dtype == bool]
    if not check_cols:
        sys.exit("Error: PoseBusters returned no boolean checks. Version mismatch?")

    df["pb_valid"] = df[check_cols].all(axis=1)

    report_path = out_stem.with_name(out_stem.name + "_report.csv")
    df.to_csv(report_path)

    passed_idx = [i for i, ok in enumerate(df["pb_valid"].tolist()) if ok]
    passed_path = out_stem.with_name(out_stem.name + "_passed.sdf")
    with Chem.SDWriter(str(passed_path)) as w:
        for i in passed_idx:
            m = mols[i]
            m.SetProp("pb_valid", "True")
            w.write(m)

    if args.keep_failed:
        failed_path = out_stem.with_name(out_stem.name + "_failed.sdf")
        with Chem.SDWriter(str(failed_path)) as w:
            for i, ok in enumerate(df["pb_valid"].tolist()):
                if not ok:
                    fails = [c for c in check_cols if not df.iloc[i][c]]
                    mols[i].SetProp("pb_failed_checks", ";".join(fails))
                    w.write(mols[i])
        print(f"failed poses  -> {failed_path}")

    n_pass = len(passed_idx)
    summary_lines = [
        f"input:  {args.poses}",
        f"mode:   {mode}",
        f"poses:  {len(mols)}",
        f"passed: {n_pass} ({100 * n_pass / len(mols):.1f}%)",
        f"failed: {len(mols) - n_pass}",
        "",
        "failures by check:",
    ]
    for col in check_cols:
        n_fail = int((~df[col]).sum())
        if n_fail:
            summary_lines.append(f"  {col:<42} {n_fail:>6}")

    summary = "\n".join(summary_lines)
    summary_path = out_stem.with_name(out_stem.name + "_summary.txt")
    summary_path.write_text(summary + "\n")

    print()
    print(summary)
    print(f"\nreport  -> {report_path}")
    print(f"passed  -> {passed_path}")

    if n_pass == 0:
        print("\nWARNING: nothing passed. Before blaming the ligands, check the "
              "receptor file and the docking box -- a misplaced box fails every "
              "pose on protein-ligand distance checks.", file=sys.stderr)

    return 0


if __name__ == "__main__":
    sys.exit(main())
