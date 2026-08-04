#!/usr/bin/env bash
set -euo pipefail

# Run only the hypothetical/uncharacterized protein analysis.
#
# Usage:
#   bash run_hypothetical.sh [protein_fasta] [Pfam-A.hmm] [cpu] [prefix]

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$ROOT"

PROTEINS="${1:-data/protein.faa}"
PFAM_DB="${2:-data/Pfam-A.hmm}"
CPU="${3:-4}"
PREFIX="${4:-ECOLI_K12}"

OUTPUT_DIR="results/hypothetical"
HYP_FASTA="${OUTPUT_DIR}/${PREFIX}_hypothetical_proteins.faa"

python3 scripts/extract_hypothetical_proteins.py \
    --input "$PROTEINS" \
    --output "$HYP_FASTA"

HYP_COUNT="$(grep -c '^>' "$HYP_FASTA" || true)"

if [[ "$HYP_COUNT" -eq 0 ]]; then
    echo "No hypothetical or uncharacterized proteins were detected."
    exit 0
fi

python3 scripts/pfam_scan.py \
    --proteins "$HYP_FASTA" \
    --pfam-db "$PFAM_DB" \
    --output-dir "$OUTPUT_DIR" \
    --prefix "${PREFIX}_hypothetical_pfama" \
    --summary-title "${PREFIX} Hypothetical-Protein Pfam-A Summary" \
    --cpu "$CPU"
