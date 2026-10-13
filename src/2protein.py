#!/usr/bin/env python3
"""
Stage 2: receptor preparation.

Reads an mmCIF or PDB structure, reports what is actually in it, strips the
things that should not be docked against (waters, cryoprotectants, buffer
ions, alternate conformers), and writes:

  <out>.pdb         cleaned receptor
  <out>_box.json    docking box, derived from the co-crystallized ligand
  <out>_ligand.sdf  the extracted native ligand (redocking control)

The box comes from the crystal ligand by default. That is the single most
important input to a docking run and the easiest thing to get silently wrong
-- a box centered on the wrong site produces confident, meaningless scores.
Use --list-het first to see what is there before trusting --autobox-ligand.

Hydrogens are NOT added here. GNINA handles protonation internally, and
adding them with the wrong tool is worse than not adding them.
"""

import argparse
import json
import sys
from pathlib import Path

import gemmi


# Residues that are structural noise rather than binding partners. Not
# exhaustive -- always eyeball --list-het output for a new structure.
JUNK = {
    "HOH", "DOD", "WAT",                      # water
    "SO4", "PO4", "NO3", "CL", "BR", "IOD",   # buffer ions
    "NA", "K", "MG", "CA", "ZN", "MN", "CD", "NI", "CO", "CU", "FE", "HG",
    "GOL", "EDO", "PEG", "PG4", "PGE", "MPD", "DMS", "TRS", "ACT", "FMT",
    "EOH", "IPA", "MES", "EPE", "BME", "DTT", "CIT", "TLA", "ACY",
}


def parse_args():
    p = argparse.ArgumentParser(description="Prepare a receptor for docking.")
    p.add_argument("input", type=Path, help="mmCIF or PDB structure")
    p.add_argument("-o", "--output", type=Path, default=None,
                   help="Output stem (default: <input stem>_receptor)")
    p.add_argument("--list-het", action="store_true",
                   help="List non-protein residues and exit. Run this first.")
    p.add_argument("--chains", default=None,
                   help="Comma-separated chains to keep (default: all protein chains)")
    p.add_argument("--autobox-ligand", default=None,
                   help="Residue name of the co-crystal ligand to build the box "
                        "around (e.g. EST). Default: the largest non-junk HET.")
    p.add_argument("--box-padding", type=float, default=4.0,
                   help="Angstroms added around the ligand extent (default: 4.0)")
    p.add_argument("--box-center", default=None,
                   help="Explicit 'x,y,z' center, overriding --autobox-ligand")
    p.add_argument("--box-size", default=None,
                   help="Explicit 'x,y,z' size in Angstroms")
    p.add_argument("--keep-het", default="",
                   help="Comma-separated HET residues to KEEP in the receptor "
                        "(cofactors, structural metals)")
    p.add_argument("--altloc", default="A",
                   help="Alternate conformer to keep (default: A)")
    return p.parse_args()


def survey(st):
    """Inventory of what is in the structure."""
    protein_chains, het = {}, {}
    for model in st:
        for chain in model:
            n_aa = 0
            for res in chain:
                info = gemmi.find_tabulated_residue(res.name)
                if info and info.is_amino_acid():
                    n_aa += 1
                else:
                    key = (res.name, chain.name)
                    het[key] = het.get(key, 0) + len(res)
            if n_aa:
                protein_chains[chain.name] = n_aa
        break  # first model only
    return protein_chains, het


def main():
    args = parse_args()

    if not args.input.is_file():
        sys.exit(f"Error: file not found: {args.input}")

    st = gemmi.read_structure(str(args.input))
    st.setup_entities()
    if args.altloc:
        st.remove_alternative_conformations()

    protein_chains, het = survey(st)

    print(f"structure: {args.input.name}")
    print(f"  name:    {st.name}")
    print(f"  spacegroup: {st.spacegroup_hm or 'n/a'}")
    print(f"  models:  {len(st)}")
    print(f"  protein chains: " +
          ", ".join(f"{c} ({n} aa)" for c, n in protein_chains.items()))

    if het:
        print("  non-protein residues:")
        for (name, chain), natoms in sorted(het.items(), key=lambda kv: -kv[1]):
            tag = "  [junk]" if name in JUNK else ""
            print(f"      {name:<6} chain {chain:<2} {natoms:>4} atoms{tag}")
    else:
        print("  non-protein residues: none")

    if args.list_het:
        return 0

    if not protein_chains:
        sys.exit("Error: no protein chains found. Is this the right file?")

    keep_het = {h.strip().upper() for h in args.keep_het.split(",") if h.strip()}
    wanted_chains = ({c.strip() for c in args.chains.split(",")}
                     if args.chains else set(protein_chains))

    # Pick the autobox ligand before stripping anything.
    ligand_name = args.autobox_ligand
    if ligand_name is None and args.box_center is None:
        candidates = {k: v for k, v in het.items()
                      if k[0] not in JUNK and k[0] not in keep_het}
        if not candidates:
            sys.exit("Error: no co-crystal ligand found for the box. Pass "
                     "--box-center/--box-size explicitly, or --autobox-ligand.")
        ligand_name = max(candidates.items(), key=lambda kv: kv[1])[0][0]
        print(f"\nautobox ligand: {ligand_name} (largest non-junk HET)")
        if len(candidates) > 1:
            print("  NOTE: more than one candidate. Verify this is the "
                  "binding site you want:")
            for (name, chain), n in sorted(candidates.items(), key=lambda kv: -kv[1]):
                print(f"      {name} (chain {chain}, {n} atoms)")

    # Extract ligand atoms for the box and for the redocking control.
    lig_atoms = []
    if ligand_name:
        for chain in st[0]:
            for res in chain:
                if res.name == ligand_name:
                    lig_atoms.extend((a.pos.x, a.pos.y, a.pos.z)
                                     for a in res if a.element.name != "H")
        if not lig_atoms:
            sys.exit(f"Error: ligand {ligand_name} not found in structure.")

    # Compute the box.
    if args.box_center:
        cx, cy, cz = (float(v) for v in args.box_center.split(","))
    else:
        cx = sum(a[0] for a in lig_atoms) / len(lig_atoms)
        cy = sum(a[1] for a in lig_atoms) / len(lig_atoms)
        cz = sum(a[2] for a in lig_atoms) / len(lig_atoms)

    if args.box_size:
        sx, sy, sz = (float(v) for v in args.box_size.split(","))
    elif lig_atoms:
        pad = args.box_padding
        sx = max(a[0] for a in lig_atoms) - min(a[0] for a in lig_atoms) + 2 * pad
        sy = max(a[1] for a in lig_atoms) - min(a[1] for a in lig_atoms) + 2 * pad
        sz = max(a[2] for a in lig_atoms) - min(a[2] for a in lig_atoms) + 2 * pad
    else:
        sys.exit("Error: need --box-size when using --box-center without a ligand.")

    out_stem = args.output or args.input.with_name(f"{args.input.stem}_receptor")
    out_stem = Path(out_stem)
    out_stem.parent.mkdir(parents=True, exist_ok=True)

    # Write the native ligand as its own file, before stripping it.
    if ligand_name:
        lig_st = gemmi.Structure()
        lig_model = gemmi.Model("1")
        lig_chain = gemmi.Chain("L")
        for chain in st[0]:
            for res in chain:
                if res.name == ligand_name:
                    lig_chain.add_residue(res)
                    break
            if len(lig_chain):
                break
        lig_model.add_chain(lig_chain)
        lig_st.add_model(lig_model)
        lig_st.setup_entities()
        lig_path = out_stem.with_name(out_stem.name + "_ligand.pdb")
        lig_st.write_pdb(str(lig_path))
        print(f"\nnative ligand -> {lig_path}")
        print("  (convert to SDF with: obabel -ipdb <in> -osdf -O <out>)")

    # Strip: unwanted chains, then non-protein residues within the kept ones.
    # gemmi's Chain.__delitem__ takes an integer index, so delete in reverse
    # to keep the remaining indices valid.
    n_removed = 0
    for model in st:
        for chain in list(model):
            if chain.name not in wanted_chains:
                model.remove_chain(chain.name)
                continue
            for i in range(len(chain) - 1, -1, -1):
                res = chain[i]
                info = gemmi.find_tabulated_residue(res.name)
                is_aa = bool(info and info.is_amino_acid())
                if not is_aa and res.name not in keep_het:
                    del chain[i]
                    n_removed += 1

    st.remove_hydrogens()
    st.remove_empty_chains()
    print(f"\nstripped {n_removed} non-protein residues"
          + (f" (kept: {', '.join(sorted(keep_het))})" if keep_het else ""))

    pdb_path = out_stem.with_suffix(".pdb")
    st.write_pdb(str(pdb_path))

    box = {
        "center_x": round(cx, 3),
        "center_y": round(cy, 3),
        "center_z": round(cz, 3),
        "size_x": round(sx, 3),
        "size_y": round(sy, 3),
        "size_z": round(sz, 3),
        "source": ("explicit" if args.box_center
                   else f"autobox:{ligand_name} pad={args.box_padding}"),
        "receptor": pdb_path.name,
        "input": args.input.name,
    }
    box_path = out_stem.with_name(out_stem.name + "_box.json")
    with open(box_path, "w") as fh:
        json.dump(box, fh, indent=2)

    print(f"receptor      -> {pdb_path}")
    print(f"box           -> {box_path}")
    print(f"\nbox center: ({cx:.2f}, {cy:.2f}, {cz:.2f})")
    print(f"box size:   ({sx:.2f}, {sy:.2f}, {sz:.2f})  "
          f"volume {sx * sy * sz:.0f} A^3")

    if max(sx, sy, sz) > 40:
        print("\nWARNING: box dimension >40 A. That is large for a single site "
              "and will slow docking and admit spurious poses.", file=sys.stderr)

    return 0


if __name__ == "__main__":
    sys.exit(main())