#!/usr/bin/env python3
"""
Identify DUF-containing proteins from a Pfam full-results CSV file.

This script expects one row per protein-Pfam pair, such as the CSV produced by
scripts/pfam_scan.py. It is also tolerant of several common alternative column
names.

DUF classification:
  - Pfam name or description contains a term such as DUF1234
  - Pfam description contains "domain of unknown function"

Outputs:
  <prefix>_all_DUF_matches.csv
  <prefix>_DUF_domain_counts.csv
  <prefix>_proteins_with_DUFs.csv
  <prefix>_DUF_summary.txt
"""

from __future__ import annotations

import argparse
import csv
import re
import sys
from collections import defaultdict
from pathlib import Path
from typing import Iterable


DUF_PATTERN = re.compile(r"(?<![A-Za-z0-9])DUF[\s_-]*\d+(?!\d)", re.IGNORECASE)
UNKNOWN_FUNCTION_PATTERN = re.compile(
    r"\bdomain(?:s)?\s+of\s+unknown\s+function\b", re.IGNORECASE
)


def fail(message: str) -> "NoReturn":
    raise SystemExit(f"ERROR: {message}")


def normalize_header(value: str) -> str:
    return re.sub(r"[^a-z0-9]+", "_", value.strip().lower()).strip("_")


def choose_column(
    normalized_to_original: dict[str, str],
    candidates: Iterable[str],
    required: bool = True,
) -> str | None:
    for candidate in candidates:
        if candidate in normalized_to_original:
            return normalized_to_original[candidate]

    if required:
        available = ", ".join(sorted(normalized_to_original.values()))
        fail(
            "Could not identify a required column. "
            f"Tried: {', '.join(candidates)}. Available columns: {available}"
        )
    return None


def is_duf(pfam_name: str, description: str) -> bool:
    combined = f"{pfam_name} {description}".strip()
    return bool(
        DUF_PATTERN.search(combined)
        or UNKNOWN_FUNCTION_PATTERN.search(combined)
    )


def write_csv(path: Path, fieldnames: list[str], rows: Iterable[dict]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Extract and summarize DUF matches from Pfam results."
    )
    parser.add_argument(
        "--input",
        required=True,
        type=Path,
        help="Full Pfam protein-Pfam results CSV.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("results/duf_analysis"),
        help="Output directory (default: results/duf_analysis).",
    )
    parser.add_argument(
        "--prefix",
        default="ECOLI_K12_pfama",
        help="Output filename prefix (default: ECOLI_K12_pfama).",
    )
    return parser


def main() -> None:
    args = build_parser().parse_args()

    input_path = args.input.expanduser().resolve()
    output_dir = args.output_dir.expanduser().resolve()

    if not input_path.is_file():
        fail(f"Input CSV not found: {input_path}")

    output_dir.mkdir(parents=True, exist_ok=True)

    with input_path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)

        if not reader.fieldnames:
            fail(f"No CSV header was detected in {input_path}")

        normalized_to_original = {
            normalize_header(name): name for name in reader.fieldnames
        }

        protein_col = choose_column(
            normalized_to_original,
            [
                "protein",
                "protein_id",
                "query",
                "query_id",
                "query_name",
                "sequence",
                "sequence_id",
            ],
        )
        accession_col = choose_column(
            normalized_to_original,
            [
                "pfam_accession",
                "pfam_acc",
                "accession",
                "hmm_accession",
                "model_accession",
            ],
        )
        name_col = choose_column(
            normalized_to_original,
            [
                "pfam_name",
                "pfamid",
                "pfam_id",
                "model",
                "model_name",
                "hmm_name",
                "target_name",
            ],
            required=False,
        )
        description_col = choose_column(
            normalized_to_original,
            [
                "description",
                "pfam_description",
                "model_description",
                "target_description",
            ],
            required=False,
        )

        original_fields = list(reader.fieldnames)
        rows = list(reader)

    # Keep one row per unique protein-Pfam pair.
    # This avoids counting repeated copies of the same DUF within one protein
    # as separate protein-DUF relationships.
    unique_pairs: dict[tuple[str, str], dict] = {}

    for row in rows:
        protein = (row.get(protein_col) or "").strip()
        accession = (row.get(accession_col) or "").strip().split(".", 1)[0]
        pfam_name = (row.get(name_col) or "").strip() if name_col else ""
        description = (
            (row.get(description_col) or "").strip()
            if description_col
            else ""
        )

        if not protein or not accession:
            continue

        if not is_duf(pfam_name, description):
            continue

        key = (protein, accession)
        if key not in unique_pairs:
            cleaned = dict(row)
            cleaned[protein_col] = protein
            cleaned[accession_col] = accession
            unique_pairs[key] = cleaned

    duf_rows = sorted(
        unique_pairs.values(),
        key=lambda row: (
            row.get(protein_col, ""),
            row.get(accession_col, ""),
        ),
    )

    proteins_by_duf: dict[str, set[str]] = defaultdict(set)
    names_by_duf: dict[str, str] = {}
    descriptions_by_duf: dict[str, str] = {}
    dufs_by_protein: dict[str, set[str]] = defaultdict(set)

    for row in duf_rows:
        protein = row[protein_col]
        accession = row[accession_col]
        pfam_name = (row.get(name_col) or "").strip() if name_col else ""
        description = (
            (row.get(description_col) or "").strip()
            if description_col
            else ""
        )

        proteins_by_duf[accession].add(protein)
        dufs_by_protein[protein].add(accession)

        if pfam_name:
            names_by_duf[accession] = pfam_name
        if description:
            descriptions_by_duf[accession] = description

    all_matches_path = output_dir / f"{args.prefix}_all_DUF_matches.csv"
    write_csv(all_matches_path, original_fields, duf_rows)

    domain_rows = []
    for accession in sorted(
        proteins_by_duf,
        key=lambda acc: (-len(proteins_by_duf[acc]), acc),
    ):
        proteins = sorted(proteins_by_duf[accession])
        domain_rows.append(
            {
                "pfam_accession": accession,
                "pfam_name": names_by_duf.get(accession, ""),
                "pfam_description": descriptions_by_duf.get(accession, ""),
                "num_proteins": len(proteins),
                "proteins": ";".join(proteins),
            }
        )

    domain_counts_path = output_dir / f"{args.prefix}_DUF_domain_counts.csv"
    write_csv(
        domain_counts_path,
        [
            "pfam_accession",
            "pfam_name",
            "pfam_description",
            "num_proteins",
            "proteins",
        ],
        domain_rows,
    )

    protein_rows = []
    for protein in sorted(
        dufs_by_protein,
        key=lambda item: (-len(dufs_by_protein[item]), item),
    ):
        accessions = sorted(dufs_by_protein[protein])
        names = [names_by_duf.get(accession, "") for accession in accessions]
        protein_rows.append(
            {
                "protein": protein,
                "num_unique_dufs": len(accessions),
                "duf_accessions": ";".join(accessions),
                "duf_names": ";".join(names),
            }
        )

    proteins_path = output_dir / f"{args.prefix}_proteins_with_DUFs.csv"
    write_csv(
        proteins_path,
        [
            "protein",
            "num_unique_dufs",
            "duf_accessions",
            "duf_names",
        ],
        protein_rows,
    )

    total_pairs = len(duf_rows)
    unique_dufs = len(proteins_by_duf)
    proteins_with_dufs = len(dufs_by_protein)
    proteins_one_duf = sum(
        1 for accessions in dufs_by_protein.values() if len(accessions) == 1
    )
    proteins_multiple_dufs = sum(
        1 for accessions in dufs_by_protein.values() if len(accessions) > 1
    )
    dufs_one_protein = sum(
        1 for proteins in proteins_by_duf.values() if len(proteins) == 1
    )
    dufs_multiple_proteins = sum(
        1 for proteins in proteins_by_duf.values() if len(proteins) > 1
    )

    summary_lines = [
        f"{args.prefix} DUF Analysis Summary",
        "================================",
        f"Input file: {input_path.name}",
        "",
        f"Total unique protein-DUF pairs: {total_pairs}",
        f"Unique DUF domains: {unique_dufs}",
        f"Proteins containing at least one DUF: {proteins_with_dufs}",
        f"Proteins containing exactly one DUF: {proteins_one_duf}",
        f"Proteins containing multiple DUFs: {proteins_multiple_dufs}",
        f"DUF domains found in exactly one protein: {dufs_one_protein}",
        f"DUF domains found in multiple proteins: {dufs_multiple_proteins}",
        "",
        "DUF definition used:",
        "- Pfam name or description contains DUF followed by a number",
        '- or description contains "domain of unknown function"',
        "",
        "Files generated:",
        f"- {all_matches_path.name}",
        f"- {domain_counts_path.name}",
        f"- {proteins_path.name}",
        f"- {args.prefix}_DUF_summary.txt",
    ]

    summary_path = output_dir / f"{args.prefix}_DUF_summary.txt"
    summary_path.write_text(
        "\n".join(summary_lines) + "\n",
        encoding="utf-8",
    )

    print("DUF analysis complete")
    print(f"  Unique protein-DUF pairs: {total_pairs}")
    print(f"  Unique DUF domains: {unique_dufs}")
    print(f"  Proteins with at least one DUF: {proteins_with_dufs}")
    print(f"  Output directory: {output_dir}")


if __name__ == "__main__":
    main()
