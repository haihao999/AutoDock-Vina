#!/usr/bin/env python3

import argparse
import csv
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import AllChem


def build_mol(smiles: str) -> Chem.Mol | None:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    mol = Chem.AddHs(mol)
    if AllChem.EmbedMolecule(mol, AllChem.ETKDG()) != 0:
        return None
    AllChem.UFFOptimizeMolecule(mol)
    return mol


def write_sdf(mol: Chem.Mol, out_path: Path) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    writer = Chem.SDWriter(str(out_path))
    writer.write(mol)
    writer.close()


def generate_sdfs(csv_path: Path, output_dir: Path, substrates: set[str]) -> list[Path]:
    written = []
    with csv_path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        for index, row in enumerate(reader, start=1):
            substrate = (row.get("Substrate") or "").strip()
            if substrates and substrate not in substrates:
                continue
            smiles = (row.get("SMILES") or "").strip()
            pdbpath = (row.get("pdbpath") or "").strip()
            if not substrate or not smiles:
                continue
            mol = build_mol(smiles)
            if mol is None:
                continue
            pdb_stem = Path(pdbpath).stem if pdbpath else f"row{index}"
            filename = f"{substrate}_{pdb_stem}.sdf"
            out_path = output_dir / filename
            write_sdf(mol, out_path)
            written.append(out_path)
    return written


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Generate 3D SDF ligands from a CSV file (Substrate/SMILES)."
    )
    parser.add_argument(
        "--input",
        required=True,
        type=Path,
        help="CSV with at least Substrate and SMILES columns.",
    )
    parser.add_argument(
        "--output-dir",
        required=True,
        type=Path,
        help="Directory to write SDF files.",
    )
    parser.add_argument(
        "--substrates",
        default="PRPn,RPnTP",
        help="Comma-separated list of substrate names to include (default: PRPn,RPnTP).",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    substrates = {s.strip() for s in args.substrates.split(",") if s.strip()}
    written = generate_sdfs(args.input, args.output_dir, substrates)
    if not written:
        print("No SDF files were generated.")
        return
    print("Generated SDF files:")
    for path in written:
        print(f"- {path}")


if __name__ == "__main__":
    main()
