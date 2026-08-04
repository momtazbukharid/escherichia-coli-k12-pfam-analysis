#!/usr/bin/env bash
set -euo pipefail

# Complete Pfam workflow for an Escherichia coli K-12 protein proteome.
#
# Usage:
#   bash run_all.sh [protein_fasta] [Pfam-A.hmm] [cpu] [output_prefix]
#
# Example:
#   bash run_all.sh data/protein.faa /path/to/Pfam-A.hmm 8 ECOLI_K12
#
# The default prefix is ECOLI_K12. For an exact derivative, use a more
# specific prefix such as ECOLI_K12_MG1655 or ECOLI_K12_W3110.

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$ROOT"

PROTEINS="${1:-data/protein.faa}"
PFAM_DB="${2:-data/Pfam-A.hmm}"
CPU="${3:-4}"
PREFIX="${4:-ECOLI_K12}"

if [[ ! -f "$PROTEINS" ]]; then
    echo "ERROR: Protein FASTA not found: $PROTEINS" >&2
    exit 1
fi

if [[ ! -f "$PFAM_DB" ]]; then
    echo "ERROR: Pfam-A database not found: $PFAM_DB" >&2
    echo "Keep Pfam-A.hmm locally; it does not need to be uploaded to GitHub." >&2
    exit 1
fi

find_list() {
    for candidate in "$@"; do
        if [[ -f "$candidate" ]]; then
            printf '%s\n' "$candidate"
            return 0
        fi
    done
    return 1
}

ESSENTIAL_LIST="$(find_list \
    lists/essential102_pfams.csv \
    lists/essential102_pfams.txt || true)"

FASTAAI_LIST="$(find_list \
    lists/fastaai122_pfams.csv \
    lists/fastaai122_pfams.txt || true)"

PERSISTENT_LIST="$(find_list \
    lists/persistent500_pfams.csv \
    lists/persistent500_pfams.txt || true)"

FULL_DIR="results/pfam-a"
FULL_RESULTS="${FULL_DIR}/${PREFIX}_pfamA_full_results.csv"

echo "============================================================"
echo "1. Complete Pfam-A scan"
echo "============================================================"

python3 scripts/pfam_scan.py \
    --proteins "$PROTEINS" \
    --pfam-db "$PFAM_DB" \
    --output-dir "$FULL_DIR" \
    --prefix "${PREFIX}_pfamA" \
    --full-results-name "${PREFIX}_pfamA_full_results.csv" \
    --summary-name "${PREFIX}_pfamA_summary.txt" \
    --summary-title "${PREFIX} Complete Pfam-A hmmscan Summary" \
    --cpu "$CPU"

if [[ -n "$ESSENTIAL_LIST" ]]; then
    echo
    echo "============================================================"
    echo "2. Essential-102 Pfam scan"
    echo "============================================================"

    python3 scripts/pfam_scan.py \
        --proteins "$PROTEINS" \
        --pfam-db "$PFAM_DB" \
        --subset-ids "$ESSENTIAL_LIST" \
        --output-dir results/102_essential \
        --prefix "${PREFIX}_essential102" \
        --summary-title "${PREFIX} Essential-102 Pfam Summary" \
        --subset-label "Essential Pfam IDs" \
        --extracted-hmm-name "${PREFIX}_Essential102_extracted.hmm" \
        --cpu "$CPU"
else
    echo "[skip] Essential-102 list not found in lists/."
fi

if [[ -n "$FASTAAI_LIST" ]]; then
    echo
    echo "============================================================"
    echo "3. FastAAI-122 Pfam scan"
    echo "============================================================"

    python3 scripts/pfam_scan.py \
        --proteins "$PROTEINS" \
        --pfam-db "$PFAM_DB" \
        --subset-ids "$FASTAAI_LIST" \
        --output-dir results/122_fastaai \
        --prefix "${PREFIX}_fastaai122" \
        --summary-title "${PREFIX} FastAAI-122 Pfam Summary" \
        --subset-label "FastAAI Pfam IDs" \
        --extracted-hmm-name "${PREFIX}_FastAAI122_extracted.hmm" \
        --cpu "$CPU"
else
    echo "[skip] FastAAI-122 list not found in lists/."
fi

if [[ -n "$PERSISTENT_LIST" ]]; then
    echo
    echo "============================================================"
    echo "4. Persistent-500 Pfam scan"
    echo "============================================================"

    python3 scripts/pfam_scan.py \
        --proteins "$PROTEINS" \
        --pfam-db "$PFAM_DB" \
        --subset-ids "$PERSISTENT_LIST" \
        --output-dir results/pfam_500 \
        --prefix "${PREFIX}_persistent500" \
        --summary-title "${PREFIX} Persistent-500 Pfam Summary" \
        --subset-label "Persistent Pfam IDs" \
        --extracted-hmm-name "${PREFIX}_Persistent500_extracted.hmm" \
        --cpu "$CPU"
else
    echo "[skip] Persistent-500 list not found in lists/."
fi

echo
echo "============================================================"
echo "5. Hypothetical/uncharacterized protein analysis"
echo "============================================================"

HYP_DIR="results/hypothetical"
HYP_FASTA="${HYP_DIR}/${PREFIX}_hypothetical_proteins.faa"

python3 scripts/extract_hypothetical_proteins.py \
    --input "$PROTEINS" \
    --output "$HYP_FASTA"

HYP_COUNT="$(grep -c '^>' "$HYP_FASTA" || true)"

if [[ "$HYP_COUNT" -gt 0 ]]; then
    python3 scripts/pfam_scan.py \
        --proteins "$HYP_FASTA" \
        --pfam-db "$PFAM_DB" \
        --output-dir "$HYP_DIR" \
        --prefix "${PREFIX}_hypothetical_pfama" \
        --summary-title "${PREFIX} Hypothetical-Protein Pfam-A Summary" \
        --cpu "$CPU"
else
    echo "[skip] No hypothetical or uncharacterized proteins were detected."
fi

echo
echo "============================================================"
echo "6. DUF analysis"
echo "============================================================"

python3 scripts/analyze_dufs.py \
    --input "$FULL_RESULTS" \
    --output-dir results/duf_analysis \
    --prefix "${PREFIX}_pfamA"

echo
echo "All available analyses finished."
echo "Results directory: $ROOT/results"
