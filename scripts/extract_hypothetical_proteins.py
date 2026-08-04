#!/usr/bin/env python3
"""Extract hypothetical/uncharacterized proteins from an annotated protein FASTA."""

from __future__ import annotations

import argparse
from pathlib import Path


def fasta_records(path: Path):
    header = None
    sequence_lines = []

    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(sequence_lines)
                header = line.rstrip("\n")
                sequence_lines = []
            else:
                sequence_lines.append(line.strip())

    if header is not None:
        yield header, "".join(sequence_lines)


def wrap(sequence: str, width: int = 60) -> str:
    return "\n".join(sequence[i:i + width] for i in range(0, len(sequence), width))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, type=Path, help="Annotated protein FASTA.")
    parser.add_argument("--output", required=True, type=Path, help="Output FASTA.")
    parser.add_argument(
        "--keywords",
        nargs="+",
        default=["hypothetical protein", "uncharacterized protein"],
        help="Case-insensitive phrases searched in FASTA headers.",
    )
    args = parser.parse_args()

    keywords = [item.lower() for item in args.keywords]
    args.output.parent.mkdir(parents=True, exist_ok=True)

    total = 0
    selected = 0
    with args.output.open("w", encoding="utf-8") as destination:
        for header, sequence in fasta_records(args.input):
            total += 1
            if any(keyword in header.lower() for keyword in keywords):
                selected += 1
                destination.write(f"{header}\n{wrap(sequence)}\n")

    print(f"Total proteins: {total}")
    print(f"Selected hypothetical/uncharacterized proteins: {selected}")
    print(f"Output: {args.output}")


if __name__ == "__main__":
    main()
