#!/usr/bin/env python3
"""Generate eCM13/SDA AlphaFold 3 inputs from Figure 5 synteny targets."""

from __future__ import annotations

import argparse
import csv
import json
import re
import textwrap
from pathlib import Path

TYPE_TO_ROLE = {
    "sdaAB (Beta)": "sdaAB_beta",
    "sdaAA (Alpha)": "sdaAA_alpha",
    "sdaA (Fused)": "sdaA_fused",
}

SPLIT_ORDER = ["sdaAB_beta", "sdaAA_alpha"]


def repo_root() -> Path:
    return Path(__file__).resolve().parents[3]


def clean_sequence(seq: str) -> str:
    return re.sub(r"\s+", "", seq).rstrip("*")


def read_fasta(path: Path) -> dict[str, str]:
    sequences: dict[str, str] = {}
    current_id: str | None = None
    chunks: list[str] = []

    with path.open("r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_id is not None:
                    sequences[current_id] = clean_sequence("".join(chunks))
                current_id = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line)

    if current_id is not None:
        sequences[current_id] = clean_sequence("".join(chunks))

    return sequences


def parse_molecule(value: str) -> tuple[str, str]:
    match = re.search(r"\((human\d+)\)", value)
    if not match:
        raise ValueError(f"Cannot parse MAG id from molecule: {value}")
    mag = match.group(1)
    species = value[: match.start()].strip()
    return species, mag


def safe_name(value: str) -> str:
    return re.sub(r"_+", "_", re.sub(r"[^A-Za-z0-9]+", "_", value)).strip("_")


def read_targets(path: Path) -> dict[str, dict[str, object]]:
    targets: dict[str, dict[str, object]] = {}
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            target_type = row["type"]
            if target_type == "blank":
                continue
            if target_type not in TYPE_TO_ROLE:
                raise ValueError(f"Unsupported target type: {target_type}")
            species, mag = parse_molecule(row["molecule"])
            role = TYPE_TO_ROLE[target_type]
            entry = targets.setdefault(mag, {"species": species, "genes": {}})
            entry["genes"][role] = row["gene"]
    return targets


def build_outputs(targets: dict[str, dict[str, object]], sequences: dict[str, str]) -> tuple[list[dict], list[tuple[str, str, str, str]]]:
    jobs: list[dict] = []
    fasta_records: list[tuple[str, str, str, str]] = []

    for mag, entry in targets.items():
        species = str(entry["species"])
        genes = dict(entry["genes"])
        roles = ["sdaA_fused"] if "sdaA_fused" in genes else SPLIT_ORDER
        missing_roles = [role for role in roles if role not in genes]
        if missing_roles:
            raise ValueError(f"Incomplete eCM13 target set for {mag}: missing {', '.join(missing_roles)}")

        chains = []
        for role in roles:
            gene = str(genes[role])
            seq_key = f"{mag}|{gene}"
            if seq_key not in sequences:
                raise ValueError(f"Missing sequence in source FASTA: {seq_key}")
            sequence = sequences[seq_key]
            if not sequence:
                raise ValueError(f"Empty sequence after cleanup: {seq_key}")
            chains.append({"proteinChain": {"sequence": sequence, "count": 1}})
            fasta_records.append((seq_key, role, species, sequence))

        architecture = "sdaA_fused" if roles == ["sdaA_fused"] else "sdaAB_sdaAA_split"
        jobs.append({
            "name": f"{safe_name(species)}_{mag}_{architecture}",
            "modelSeeds": [],
            "sequences": chains,
            "dialect": "alphafoldserver",
            "version": 1,
        })

    return jobs, fasta_records


def write_fasta(records: list[tuple[str, str, str, str]], path: Path) -> None:
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        for seq_id, role, species, sequence in records:
            handle.write(f">{seq_id} {role} {species}\n")
            handle.write("\n".join(textwrap.wrap(sequence, width=80)))
            handle.write("\n")


def main() -> None:
    root = repo_root()
    default_synteny = root / "data" / "figure5" / "sda_synteny_clean_core_only_scaled.tsv"
    default_output_dir = Path(__file__).resolve().parent
    default_source = root / "data" / "external" / "m13_targets_combined.faa"

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--synteny", type=Path, default=default_synteny)
    parser.add_argument("--source-fasta", type=Path, default=default_source)
    parser.add_argument("--output-dir", type=Path, default=default_output_dir)
    args = parser.parse_args()

    targets = read_targets(args.synteny)
    sequences = read_fasta(args.source_fasta)
    jobs, fasta_records = build_outputs(targets, sequences)

    args.output_dir.mkdir(parents=True, exist_ok=True)
    fasta_path = args.output_dir / "eCM13_sda_sequences.fasta"
    json_path = args.output_dir / "eCM13_sda_AF3_batch.json"

    write_fasta(fasta_records, fasta_path)
    with json_path.open("w", encoding="utf-8", newline="\n") as handle:
        json.dump(jobs, handle, indent=2)
        handle.write("\n")

    print(f"Wrote {len(fasta_records)} FASTA records to {fasta_path}")
    print(f"Wrote {len(jobs)} AF3 jobs to {json_path}")


if __name__ == "__main__":
    main()
