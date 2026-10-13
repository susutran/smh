#!/usr/bin/env python3
"""
Stage 0: cheap 2D property filters.

Runs BEFORE conformer generation. Every molecule removed here is one you do
not pay ETKDG embedding time for, and 2D descriptors are ~4 orders of
magnitude cheaper than 3D embedding.

Reads a whitespace- or tab-delimited SMILES file, applies configurable
property ranges and structural alerts, writes a filtered file in the same
format plus a per-rule rejection report.

Thresholds are configurable via --config (JSON) so the screen's criteria are
recorded as data, not buried in code. Defaults are lead-like, deliberately
looser than Lipinski on MW because docking hit rates improve with some
headroom.
"""

import argparse
import json
import os
import sys
import time
from collections import Counter
from pathlib import Path

from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors, Crippen, rdMolDescriptors
from rdkit.Chem.FilterCatalog import FilterCatalog, FilterCatalogParams
from rdkit.Chem.MolStandardize import rdMolStandardize

RDLogger.DisableLog("rdApp.*")


# Lead-like defaults. Override with --config.
DEFAULTS = {
    "mw_min": 150.0,
    "mw_max": 500.0,
    "logp_min": -1.0,
    "logp_max": 5.0,
    "hbd_max": 5,
    "hba_max": 10,
    "tpsa_min": 0.0,
    "tpsa_max": 140.0,
    "rotb_max": 10,
    "rings_min": 1,
    "heavy_atoms_min": 10,
    "heavy_atoms_max": 50,
    "formal_charge_abs_max": 2,
    "allowed_elements": ["H", "C", "N", "O", "F", "P", "S", "Cl", "Br", "I"],
    "pains": True,
    "brenk": False,
    "nih": False,
    "strip_salts": True,
    "require_single_fragment": True,
}


def parse_args():
    p = argparse.ArgumentParser(description="Cheap 2D filters for virtual screening.")
    p.add_argument("input", type=Path, help="Whitespace-delimited SMILES file")
    p.add_argument("-o", "--output", type=Path, default=None,
                   help="Filtered SMILES file (default: <stem>_filtered.txt)")
    p.add_argument("--config", type=Path, default=None,
                   help="JSON file overriding any of the default thresholds")
    p.add_argument("--smiles-col", default="SMILES")
    p.add_argument("--id-col", default="zincid")
    p.add_argument("--dump-config", action="store_true",
                   help="Print the effective config as JSON and exit")
    p.add_argument("--chunk-size", type=int,
                   default=int(os.environ.get("SMH_CHUNK_SIZE", 0)))
    p.add_argument("--task-id", type=int,
                   default=int(os.environ.get("SLURM_ARRAY_TASK_ID", -1)))
    return p.parse_args()


def load_config(path):
    cfg = dict(DEFAULTS)
    if path:
        with open(path) as fh:
            user = json.load(fh)
        unknown = set(user) - set(DEFAULTS)
        if unknown:
            sys.exit(f"Error: unknown config keys: {sorted(unknown)}")
        cfg.update(user)
    return cfg


def build_alert_catalog(cfg):
    """Structural-alert catalogs. Returns None if none are enabled."""
    params = FilterCatalogParams()
    enabled = []
    if cfg["pains"]:
        params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS)
        enabled.append("PAINS")
    if cfg["brenk"]:
        params.AddCatalog(FilterCatalogParams.FilterCatalogs.BRENK)
        enabled.append("BRENK")
    if cfg["nih"]:
        params.AddCatalog(FilterCatalogParams.FilterCatalogs.NIH)
        enabled.append("NIH")
    return (FilterCatalog(params), enabled) if enabled else (None, [])


def read_rows(path, smiles_col, id_col):
    with open(path) as fh:
        header_line = fh.readline()
        delim = "\t" if "\t" in header_line else None
        header = header_line.split(delim) if delim else header_line.split()
        header = [h.strip() for h in header]
        try:
            si, ii = header.index(smiles_col), header.index(id_col)
        except ValueError:
            sys.exit(f"Error: columns {smiles_col!r}/{id_col!r} not in header: {header}")

        for lineno, line in enumerate(fh, start=2):
            parts = line.split(delim) if delim else line.split()
            parts = [p.strip() for p in parts]
            if len(parts) <= max(si, ii):
                continue
            yield lineno, parts, parts[ii], parts[si]
    return


def check(mol, cfg, catalog, allowed):
    """Return None if the molecule passes, else the name of the first failed rule."""
    if cfg["require_single_fragment"] and "." in Chem.MolToSmiles(mol):
        return "multi_fragment"

    for atom in mol.GetAtoms():
        if atom.GetSymbol() not in allowed:
            return "element"

    nheavy = mol.GetNumHeavyAtoms()
    if nheavy < cfg["heavy_atoms_min"]:
        return "heavy_atoms_min"
    if nheavy > cfg["heavy_atoms_max"]:
        return "heavy_atoms_max"

    mw = Descriptors.MolWt(mol)
    if mw < cfg["mw_min"]:
        return "mw_min"
    if mw > cfg["mw_max"]:
        return "mw_max"

    logp = Crippen.MolLogP(mol)
    if logp < cfg["logp_min"]:
        return "logp_min"
    if logp > cfg["logp_max"]:
        return "logp_max"

    if rdMolDescriptors.CalcNumHBD(mol) > cfg["hbd_max"]:
        return "hbd"
    if rdMolDescriptors.CalcNumHBA(mol) > cfg["hba_max"]:
        return "hba"

    tpsa = rdMolDescriptors.CalcTPSA(mol)
    if tpsa < cfg["tpsa_min"]:
        return "tpsa_min"
    if tpsa > cfg["tpsa_max"]:
        return "tpsa_max"

    if rdMolDescriptors.CalcNumRotatableBonds(mol) > cfg["rotb_max"]:
        return "rotatable_bonds"
    if rdMolDescriptors.CalcNumRings(mol) < cfg["rings_min"]:
        return "rings_min"

    if abs(Chem.GetFormalCharge(mol)) > cfg["formal_charge_abs_max"]:
        return "formal_charge"

    if catalog is not None and catalog.HasMatch(mol):
        entry = catalog.GetFirstMatch(mol)
        return f"alert:{entry.GetDescription()}" if entry else "alert"

    return None


def main():
    args = parse_args()
    cfg = load_config(args.config)

    if args.dump_config:
        print(json.dumps(cfg, indent=2))
        return 0

    if not args.input.is_file():
        sys.exit(f"Error: file not found: {args.input}")

    catalog, alerts_enabled = build_alert_catalog(cfg)
    allowed = set(cfg["allowed_elements"])
    remover = rdMolStandardize.LargestFragmentChooser() if cfg["strip_salts"] else None

    rows = list(read_rows(args.input, args.smiles_col, args.id_col))
    total = len(rows)

    suffix = ""
    if args.chunk_size > 0 and args.task_id >= 0:
        start = args.task_id * args.chunk_size
        rows = rows[start:start + args.chunk_size]
        suffix = f"_{args.task_id:05d}"
        if not rows:
            print(f"Task {args.task_id}: nothing in range. Done.")
            return 0

    out_path = args.output or args.input.with_name(
        f"{args.input.stem}_filtered{suffix}.txt")
    out_path.parent.mkdir(parents=True, exist_ok=True)
    rej_path = out_path.with_name(out_path.stem + ".rejected.txt")
    cfg_path = out_path.with_name(out_path.stem + ".config.json")

    # Record the exact criteria used, next to the output. Reproducibility:
    # the filtered set is meaningless without the rules that produced it.
    with open(cfg_path, "w") as fh:
        json.dump(cfg, fh, indent=2)

    print(f"input:   {args.input} ({total} rows)")
    print(f"this run:{len(rows)} molecules")
    print(f"output:  {out_path}")
    print(f"alerts:  {', '.join(alerts_enabled) if alerts_enabled else 'none'}")
    print(f"config:  {cfg_path}", flush=True)

    reasons = Counter()
    n_pass = 0
    t0 = time.time()

    with open(out_path, "w") as out, open(rej_path, "w") as rej:
        out.write(f"{args.id_col}\t{args.smiles_col}\tMW\tcLogP\tTPSA\n")
        rej.write(f"{args.id_col}\t{args.smiles_col}\treason\n")

        for i, (lineno, parts, mol_id, smiles) in enumerate(rows, start=1):
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                reasons["parse"] += 1
                rej.write(f"{mol_id}\t{smiles}\tparse\n")
                continue

            if remover is not None:
                mol = remover.choose(mol)

            reason = check(mol, cfg, catalog, allowed)
            if reason:
                reasons[reason] += 1
                rej.write(f"{mol_id}\t{smiles}\t{reason}\n")
                continue

            canonical = Chem.MolToSmiles(mol)
            out.write(f"{mol_id}\t{canonical}\t"
                      f"{Descriptors.MolWt(mol):.1f}\t"
                      f"{Crippen.MolLogP(mol):.2f}\t"
                      f"{rdMolDescriptors.CalcTPSA(mol):.1f}\n")
            n_pass += 1

            if i % 10000 == 0:
                print(f"  {i}/{len(rows)}  ({i / (time.time() - t0):.0f} mol/s)",
                      flush=True)

    elapsed = max(time.time() - t0, 1e-9)
    print(f"\ndone in {elapsed:.1f}s ({len(rows) / elapsed:.0f} mol/s)")
    print(f"  passed:   {n_pass} ({100 * n_pass / max(len(rows), 1):.1f}%)")
    print(f"  rejected: {sum(reasons.values())}")
    for reason, n in reasons.most_common():
        print(f"      {reason:<28} {n}")

    # A pass rate near 0 or 100 almost always means misconfiguration, not a
    # remarkable library. Say so rather than letting it pass silently.
    rate = n_pass / max(len(rows), 1)
    if rate < 0.01:
        print("\nWARNING: <1% passed. Check thresholds against your library.",
              file=sys.stderr)
    elif rate > 0.99 and len(rows) > 100:
        print("\nWARNING: >99% passed. Filters may not be doing anything.",
              file=sys.stderr)

    return 0 if n_pass else 1


if __name__ == "__main__":
    sys.exit(main())
