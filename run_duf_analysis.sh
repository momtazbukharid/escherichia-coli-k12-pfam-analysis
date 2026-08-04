#!/usr/bin/env bash
set -euo pipefail

# Run only the DUF analysis.
#
# Usage:
#   bash run_duf_analysis.sh [full_results.csv] [output_dir] [prefix]
#
# Example:
#   bash run_duf_analysis.sh \
#     results/pfam-a/ECOLI_K12_pfamA_full_results.csv \
#     results/duf_analysis \
#     ECOLI_K12_pfamA

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$ROOT"

INPUT="${1:-}"
OUTPUT_DIR="${2:-results/duf_analysis}"
PREFIX="${3:-ECOLI_K12_pfamA}"

if [[ -z "$INPUT" ]]; then
    CANDIDATES=(
        "results/pfam-a/${PREFIX}_full_results.csv"
        "results/pfam-a/ECOLI_K12_pfamA_full_results.csv"
        "results/pfam-a/pfamA_full_results.csv"
        "pfamA_full_results.csv"
    )

    for candidate in "${CANDIDATES[@]}"; do
        if [[ -f "$candidate" ]]; then
            INPUT="$candidate"
            break
        fi
    done
fi

if [[ -z "$INPUT" || ! -f "$INPUT" ]]; then
    echo "ERROR: Complete Pfam result CSV not found." >&2
    echo "Usage: bash run_duf_analysis.sh FULL_RESULTS.csv [OUTPUT_DIR] [PREFIX]" >&2
    exit 1
fi

python3 scripts/analyze_dufs.py \
    --input "$INPUT" \
    --output-dir "$OUTPUT_DIR" \
    --prefix "$PREFIX"
