#!/usr/bin/env python3
"""
Run HMMER hmmscan against the complete Pfam-A database or a selected Pfam subset.

The program:
  1. optionally extracts requested Pfam models from Pfam-A.hmm;
  2. presses the HMM database with hmmpress;
  3. runs hmmscan;
  4. collapses repeated domain instances to unique protein-Pfam pairs;
  5. writes full, single/multiple-Pfam, no-match, and summary reports.

Only the Python standard library is required. HMMER must be installed and
available on PATH (hmmscan and hmmpress).
"""

from __future__ import annotations

import argparse
import csv
import re
import shutil
import subprocess
import sys
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable


PFAM_RE = re.compile(r"PF\d+", re.IGNORECASE)


@dataclass
class DomainHit:
    protein: str
    pfam_accession: str
    pfam_name: str
    full_evalue: float
    full_score: float
    domain_cevalue: float
    domain_ievalue: float
    domain_score: float
    hmm_from: int
    hmm_to: int
    ali_from: int
    ali_to: int
    env_from: int
    env_to: int
    accuracy: float
    description: str
    positions: list[str] = field(default_factory=list)
    domain_instances: int = 1


def fail(message: str) -> "NoReturn":
    raise SystemExit(f"ERROR: {message}")


def require_program(name: str) -> None:
    if shutil.which(name) is None:
        fail(
            f"Required program '{name}' was not found on PATH. "
            "Install HMMER first."
        )


def run_command(command: list[str]) -> None:
    printable = " ".join(str(item) for item in command)
    print(f"\n[run] {printable}", flush=True)
    try:
        subprocess.run(command, check=True)
    except subprocess.CalledProcessError as exc:
        fail(f"Command failed with exit code {exc.returncode}: {printable}")


def normalize_pfam_id(value: str) -> str:
    """Remove an optional Pfam release suffix, e.g. PF00001.27 -> PF00001."""
    return value.strip().upper().split(".", 1)[0]


def read_fasta_ids(fasta_path: Path) -> list[str]:
    ids: list[str] = []
    seen: set[str] = set()

    with fasta_path.open("r", encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.startswith(">"):
                continue
            identifier = line[1:].strip().split()[0]
            if not identifier:
                fail(f"Empty FASTA identifier at line {line_number}.")
            if identifier in seen:
                fail(f"Duplicate FASTA identifier: {identifier}")
            seen.add(identifier)
            ids.append(identifier)

    if not ids:
        fail(f"No FASTA records were detected in {fasta_path}")
    return ids


def read_subset_ids(path: Path) -> list[str]:
    """
    Read Pfam IDs from a text/CSV/TSV file.

    Any token matching PF followed by digits is accepted. This deliberately
    preserves malformed IDs such as PF0054 so that they are reported as
    missing rather than silently discarded.
    """
    ordered: list[str] = []
    seen: set[str] = set()

    with path.open("r", encoding="utf-8-sig") as handle:
        for line in handle:
            if line.lstrip().startswith("#"):
                continue
            for match in PFAM_RE.findall(line):
                pfam_id = normalize_pfam_id(match)
                if pfam_id not in seen:
                    seen.add(pfam_id)
                    ordered.append(pfam_id)

    if not ordered:
        fail(f"No Pfam accessions were detected in subset file: {path}")
    return ordered


def extract_subset_hmms(
    pfam_db: Path,
    requested_ids: list[str],
    output_hmm: Path,
) -> tuple[list[str], list[str]]:
    """Stream through Pfam-A.hmm and copy matching HMM records."""
    requested = set(requested_ids)
    found: set[str] = set()
    output_hmm.parent.mkdir(parents=True, exist_ok=True)

    current_lines: list[str] = []
    current_accession: str | None = None

    def flush_record(out_handle) -> None:
        nonlocal current_lines, current_accession
        if current_lines and current_accession in requested:
            out_handle.writelines(current_lines)
            found.add(current_accession)
        current_lines = []
        current_accession = None

    with pfam_db.open("r", encoding="utf-8", errors="replace") as source, \
         output_hmm.open("w", encoding="utf-8") as destination:
        for line in source:
            current_lines.append(line)

            if line.startswith("ACC"):
                fields = line.split()
                if len(fields) >= 2:
                    current_accession = normalize_pfam_id(fields[1])

            if line.strip() == "//":
                flush_record(destination)

        # Defensive handling for an incomplete final record.
        if current_lines:
            flush_record(destination)

    found_ordered = [pfam_id for pfam_id in requested_ids if pfam_id in found]
    missing_ordered = [pfam_id for pfam_id in requested_ids if pfam_id not in found]

    if not found_ordered:
        fail(
            "None of the requested Pfam models were found in the supplied "
            f"database: {pfam_db}"
        )

    return found_ordered, missing_ordered


def pressed_files(hmm_path: Path) -> list[Path]:
    return [Path(str(hmm_path) + suffix) for suffix in (".h3f", ".h3i", ".h3m", ".h3p")]


def ensure_pressed(hmm_path: Path, force: bool = False) -> None:
    indices = pressed_files(hmm_path)
    if not force and all(path.exists() for path in indices):
        print(f"[skip] Existing hmmpress indexes found for {hmm_path}")
        return

    for path in indices:
        if path.exists():
            path.unlink()

    run_command(["hmmpress", "-f", str(hmm_path)])


def parse_domtblout(path: Path) -> list[DomainHit]:
    """
    Parse hmmscan --domtblout output.

    hmmscan orientation:
      target = HMM/Pfam model
      query  = protein sequence
    """
    domain_rows: list[DomainHit] = []

    with path.open("r", encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip() or line.startswith("#"):
                continue

            fields = line.rstrip("\n").split(maxsplit=22)
            if len(fields) < 22:
                fail(
                    f"Unexpected domtblout structure at line {line_number} "
                    f"in {path}"
                )

            pfam_name = fields[0]
            raw_accession = fields[1]
            protein = fields[3]
            description = fields[22] if len(fields) > 22 else ""

            if raw_accession == "-":
                pfam_accession = pfam_name
            else:
                pfam_accession = normalize_pfam_id(raw_accession)

            domain_rows.append(
                DomainHit(
                    protein=protein,
                    pfam_accession=pfam_accession,
                    pfam_name=pfam_name,
                    full_evalue=float(fields[6]),
                    full_score=float(fields[7]),
                    domain_cevalue=float(fields[11]),
                    domain_ievalue=float(fields[12]),
                    domain_score=float(fields[13]),
                    hmm_from=int(fields[15]),
                    hmm_to=int(fields[16]),
                    ali_from=int(fields[17]),
                    ali_to=int(fields[18]),
                    env_from=int(fields[19]),
                    env_to=int(fields[20]),
                    accuracy=float(fields[21]),
                    description=description,
                    positions=[f"{fields[17]}-{fields[18]}"],
                )
            )

    return domain_rows


def collapse_to_unique_pairs(domain_hits: Iterable[DomainHit]) -> list[DomainHit]:
    """
    Collapse multiple domain instances to one protein-Pfam row.

    The best instance is selected by the lowest independent E-value and then
    the highest domain score. All alignment coordinates are retained.
    """
    grouped: dict[tuple[str, str], list[DomainHit]] = defaultdict(list)
    for hit in domain_hits:
        grouped[(hit.protein, hit.pfam_accession)].append(hit)

    collapsed: list[DomainHit] = []
    for hits in grouped.values():
        best = min(hits, key=lambda h: (h.domain_ievalue, -h.domain_score))
        best.domain_instances = len(hits)
        best.positions = [f"{h.ali_from}-{h.ali_to}" for h in hits]
        collapsed.append(best)

    return collapsed


def prefixed(prefix: str, suffix: str) -> str:
    return f"{prefix}_{suffix}" if prefix else suffix


def write_csv(path: Path, fieldnames: list[str], rows: Iterable[dict]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def analyze_and_write(
    *,
    fasta_ids: list[str],
    pair_hits: list[DomainHit],
    output_dir: Path,
    prefix: str,
    bare_category_files: bool,
    summary_title: str,
    subset_label: str | None,
    requested_ids: list[str] | None,
    found_ids: list[str] | None,
    missing_ids: list[str] | None,
    extracted_hmm_name: str | None,
    tblout_name: str,
    domtblout_name: str,
    alignments_name: str,
    full_results_filename: str | None,
    summary_filename: str | None,
) -> dict[str, int]:
    output_dir.mkdir(parents=True, exist_ok=True)

    fasta_order = {protein: index for index, protein in enumerate(fasta_ids)}
    pair_hits.sort(
        key=lambda hit: (
            fasta_order.get(hit.protein, len(fasta_order)),
            hit.pfam_accession,
        )
    )

    hits_by_protein: dict[str, list[DomainHit]] = defaultdict(list)
    proteins_by_pfam: dict[str, set[str]] = defaultdict(set)
    pfam_names: dict[str, str] = {}

    for hit in pair_hits:
        hits_by_protein[hit.protein].append(hit)
        proteins_by_pfam[hit.pfam_accession].add(hit.protein)
        pfam_names[hit.pfam_accession] = hit.pfam_name

    matched_proteins = [protein for protein in fasta_ids if protein in hits_by_protein]
    unmatched_proteins = [protein for protein in fasta_ids if protein not in hits_by_protein]

    single_proteins = [
        protein for protein in matched_proteins if len(hits_by_protein[protein]) == 1
    ]
    multiple_proteins = [
        protein for protein in matched_proteins if len(hits_by_protein[protein]) > 1
    ]

    single_pfams = sorted(
        pfam for pfam, proteins in proteins_by_pfam.items() if len(proteins) == 1
    )
    multiple_pfams = sorted(
        pfam for pfam, proteins in proteins_by_pfam.items() if len(proteins) > 1
    )

    category_prefix = "" if bare_category_files else prefix

    full_results_name = full_results_filename or prefixed(prefix, "pfam_full_results.csv")
    single_proteins_name = prefixed(category_prefix, "proteins_single_pfam.csv")
    multiple_proteins_name = prefixed(category_prefix, "proteins_multiple_pfams.csv")
    single_pfams_name = prefixed(category_prefix, "pfams_single_protein.csv")
    multiple_pfams_name = prefixed(category_prefix, "pfams_multiple_proteins.csv")
    no_match_name = prefixed(category_prefix, "proteins_no_pfam_match.csv")
    summary_name = summary_filename or prefixed(prefix, "summary_report.txt")

    full_fields = [
        "protein",
        "pfam_accession",
        "pfam_name",
        "domain_instances",
        "alignment_positions",
        "full_evalue",
        "full_score",
        "best_domain_cevalue",
        "best_domain_ievalue",
        "best_domain_score",
        "hmm_from",
        "hmm_to",
        "ali_from",
        "ali_to",
        "env_from",
        "env_to",
        "accuracy",
        "description",
    ]

    full_rows = []
    for hit in pair_hits:
        full_rows.append(
            {
                "protein": hit.protein,
                "pfam_accession": hit.pfam_accession,
                "pfam_name": hit.pfam_name,
                "domain_instances": hit.domain_instances,
                "alignment_positions": ";".join(hit.positions),
                "full_evalue": hit.full_evalue,
                "full_score": hit.full_score,
                "best_domain_cevalue": hit.domain_cevalue,
                "best_domain_ievalue": hit.domain_ievalue,
                "best_domain_score": hit.domain_score,
                "hmm_from": hit.hmm_from,
                "hmm_to": hit.hmm_to,
                "ali_from": hit.ali_from,
                "ali_to": hit.ali_to,
                "env_from": hit.env_from,
                "env_to": hit.env_to,
                "accuracy": hit.accuracy,
                "description": hit.description,
            }
        )
    write_csv(output_dir / full_results_name, full_fields, full_rows)

    protein_category_fields = ["protein", "num_pfams", "pfam_accession", "pfam_name"]

    def protein_rows(proteins: list[str]) -> Iterable[dict]:
        for protein in proteins:
            hits = sorted(hits_by_protein[protein], key=lambda item: item.pfam_accession)
            for hit in hits:
                yield {
                    "protein": protein,
                    "num_pfams": len(hits),
                    "pfam_accession": hit.pfam_accession,
                    "pfam_name": hit.pfam_name,
                }

    write_csv(
        output_dir / single_proteins_name,
        protein_category_fields,
        protein_rows(single_proteins),
    )
    write_csv(
        output_dir / multiple_proteins_name,
        protein_category_fields,
        protein_rows(multiple_proteins),
    )

    pfam_category_fields = ["pfam_accession", "pfam_name", "num_proteins"]

    def pfam_rows(pfams: list[str]) -> Iterable[dict]:
        for pfam in pfams:
            yield {
                "pfam_accession": pfam,
                "pfam_name": pfam_names[pfam],
                "num_proteins": len(proteins_by_pfam[pfam]),
            }

    write_csv(output_dir / single_pfams_name, pfam_category_fields, pfam_rows(single_pfams))
    write_csv(
        output_dir / multiple_pfams_name,
        pfam_category_fields,
        pfam_rows(multiple_pfams),
    )
    write_csv(
        output_dir / no_match_name,
        ["protein"],
        ({"protein": protein} for protein in unmatched_proteins),
    )

    statistics = {
        "total_proteins": len(fasta_ids),
        "matched_proteins": len(matched_proteins),
        "unmatched_proteins": len(unmatched_proteins),
        "unique_pairs": len(pair_hits),
        "unique_pfams": len(proteins_by_pfam),
        "single_pfam_proteins": len(single_proteins),
        "multiple_pfam_proteins": len(multiple_proteins),
        "single_protein_pfams": len(single_pfams),
        "multiple_protein_pfams": len(multiple_pfams),
    }

    lines = [summary_title, "=" * len(summary_title)]

    if requested_ids is None:
        lines.extend(
            [
                f"Total proteins: {statistics['total_proteins']}",
                "",
                f"Proteins with match: {statistics['matched_proteins']}",
                f"Proteins without match: {statistics['unmatched_proteins']}",
                "",
                f"Total protein-Pfam pairs: {statistics['unique_pairs']}",
                f"Unique Pfam domains detected: {statistics['unique_pfams']}",
                "",
                f"Proteins with exactly one Pfam: {statistics['single_pfam_proteins']}",
                f"Proteins with multiple Pfams: {statistics['multiple_pfam_proteins']}",
                "",
                f"Pfams with single protein: {statistics['single_protein_pfams']}",
                f"Pfams with multiple proteins: {statistics['multiple_protein_pfams']}",
            ]
        )
    else:
        label = subset_label or "Pfam IDs"
        lines.extend(
            [
                f"Total proteins in FASTA: {statistics['total_proteins']}",
                f"{label} requested from subset file: {len(requested_ids)}",
                f"{label} successfully extracted from Pfam-A.hmm: {len(found_ids or [])}",
                f"{label} missing from Pfam-A.hmm: {len(missing_ids or [])}",
                f"Proteins with at least one Pfam match: {statistics['matched_proteins']}",
                f"Proteins without any Pfam match: {statistics['unmatched_proteins']}",
                f"Total unique protein-Pfam pairs: {statistics['unique_pairs']}",
                f"Unique Pfam domains matched: {statistics['unique_pfams']}",
                f"Proteins matching exactly one Pfam: {statistics['single_pfam_proteins']}",
                f"Proteins matching multiple Pfams: {statistics['multiple_pfam_proteins']}",
                f"Pfam domains matching exactly one protein: {statistics['single_protein_pfams']}",
                f"Pfam domains matching multiple proteins: {statistics['multiple_protein_pfams']}",
                "",
                "Missing requested Pfams from Pfam-A.hmm:",
                *(missing_ids or ["None"]),
                "",
                "Files generated:",
                f"- {full_results_name}",
                f"- {single_proteins_name}",
                f"- {multiple_proteins_name}",
                f"- {multiple_pfams_name}",
                f"- {single_pfams_name}",
                f"- {no_match_name}",
                f"- {summary_name}",
                *([f"- {extracted_hmm_name}"] if extracted_hmm_name else []),
                f"- {tblout_name}",
                f"- {domtblout_name}",
                f"- {alignments_name}",
            ]
        )

    (output_dir / summary_name).write_text("\n".join(lines) + "\n", encoding="utf-8")
    return statistics


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Run and summarize a Pfam hmmscan analysis."
    )
    parser.add_argument("--proteins", required=True, type=Path, help="Protein FASTA file.")
    parser.add_argument(
        "--pfam-db",
        required=True,
        type=Path,
        help="Complete Pfam-A.hmm database.",
    )
    parser.add_argument(
        "--subset-ids",
        type=Path,
        help="Optional file containing selected Pfam accessions.",
    )
    parser.add_argument(
        "--output-dir",
        required=True,
        type=Path,
        help="Directory for results.",
    )
    parser.add_argument(
        "--prefix",
        required=True,
        help="Prefix used for result filenames, e.g. essential102.",
    )
    parser.add_argument(
        "--summary-title",
        required=True,
        help="Title written at the top of the summary report.",
    )
    parser.add_argument(
        "--subset-label",
        help="Wording used in subset reports, e.g. 'Essential Pfam IDs'.",
    )
    parser.add_argument(
        "--extracted-hmm-name",
        help="Filename for an extracted subset HMM database.",
    )
    parser.add_argument("--cpu", type=int, default=4, help="hmmscan CPU threads.")
    parser.add_argument(
        "--threshold",
        choices=("gathering", "evalue"),
        default="gathering",
        help="Use Pfam gathering thresholds or a user-supplied E-value.",
    )
    parser.add_argument(
        "--evalue",
        type=float,
        default=1e-5,
        help="Sequence and domain E-value when --threshold evalue is selected.",
    )
    parser.add_argument(
        "--full-results-name",
        help="Optional exact filename for the full protein-Pfam results CSV.",
    )
    parser.add_argument(
        "--summary-name",
        help="Optional exact filename for the summary report.",
    )
    parser.add_argument(
        "--bare-category-files",
        action="store_true",
        help="Do not prefix category CSV filenames (matches the old full-Pfam run).",
    )
    parser.add_argument(
        "--force-press",
        action="store_true",
        help="Rebuild hmmpress indexes even when all four index files exist.",
    )
    return parser


def main() -> None:
    args = build_parser().parse_args()

    if args.cpu < 1:
        fail("--cpu must be at least 1")

    args.proteins = args.proteins.expanduser().resolve()
    args.pfam_db = args.pfam_db.expanduser().resolve()
    args.output_dir = args.output_dir.expanduser().resolve()
    if args.subset_ids:
        args.subset_ids = args.subset_ids.expanduser().resolve()

    if not args.proteins.is_file():
        fail(f"Protein FASTA not found: {args.proteins}")
    if not args.pfam_db.is_file():
        fail(f"Pfam database not found: {args.pfam_db}")
    if args.subset_ids and not args.subset_ids.is_file():
        fail(f"Subset-ID file not found: {args.subset_ids}")

    require_program("hmmscan")
    require_program("hmmpress")

    args.output_dir.mkdir(parents=True, exist_ok=True)
    fasta_ids = read_fasta_ids(args.proteins)

    requested_ids: list[str] | None = None
    found_ids: list[str] | None = None
    missing_ids: list[str] | None = None
    extracted_hmm_name: str | None = None

    if args.subset_ids:
        requested_ids = read_subset_ids(args.subset_ids)
        extracted_hmm_name = (
            args.extracted_hmm_name or f"{args.prefix}_extracted.hmm"
        )
        scan_hmm = args.output_dir / extracted_hmm_name
        found_ids, missing_ids = extract_subset_hmms(
            args.pfam_db,
            requested_ids,
            scan_hmm,
        )
    else:
        scan_hmm = args.pfam_db

    ensure_pressed(scan_hmm, force=args.force_press)

    tblout_name = f"{args.prefix}_hmmscan_tblout.txt"
    domtblout_name = f"{args.prefix}_hmmscan_domtblout.txt"
    alignments_name = f"{args.prefix}_hmmscan_alignments.txt"

    command = [
        "hmmscan",
        "--cpu",
        str(args.cpu),
        "--tblout",
        str(args.output_dir / tblout_name),
        "--domtblout",
        str(args.output_dir / domtblout_name),
        "-o",
        str(args.output_dir / alignments_name),
    ]

    if args.threshold == "gathering":
        command.append("--cut_ga")
    else:
        command.extend(["-E", str(args.evalue), "--domE", str(args.evalue)])

    command.extend([str(scan_hmm), str(args.proteins)])
    run_command(command)

    domain_hits = parse_domtblout(args.output_dir / domtblout_name)
    pair_hits = collapse_to_unique_pairs(domain_hits)

    statistics = analyze_and_write(
        fasta_ids=fasta_ids,
        pair_hits=pair_hits,
        output_dir=args.output_dir,
        prefix=args.prefix,
        bare_category_files=args.bare_category_files,
        summary_title=args.summary_title,
        subset_label=args.subset_label,
        requested_ids=requested_ids,
        found_ids=found_ids,
        missing_ids=missing_ids,
        extracted_hmm_name=extracted_hmm_name,
        tblout_name=tblout_name,
        domtblout_name=domtblout_name,
        alignments_name=alignments_name,
        full_results_filename=args.full_results_name,
        summary_filename=args.summary_name,
    )

    print("\nAnalysis complete")
    for key, value in statistics.items():
        print(f"  {key}: {value}")


if __name__ == "__main__":
    main()
