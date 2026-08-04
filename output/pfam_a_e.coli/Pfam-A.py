import subprocess
import pandas as pd
from pathlib import Path

# =========================
# INPUT FILES
# =========================
HMM_FILE = "Pfam-A.hmm"
FASTA_FILE = "protein.faa"

# =========================
# OUTPUT FILES
# =========================
TBL_OUT = "pfamA_full_tblout.txt"
ALIGNMENT_OUT = "pfamA_full_alignments.txt"
PFAM_FULL_RESULTS = "pfamA_full_results.csv"
PROTEINS_SINGLE = "pfamA_proteins_single_pfam.csv"
PROTEINS_MULTI = "pfamA_proteins_multiple_pfams.csv"
PFAMS_MULTI_PROTEINS = "pfamA_pfams_multiple_proteins.csv"
PFAMS_SINGLE_PROTEIN = "pfamA_pfams_single_protein.csv"
PROTEINS_NO_MATCH = "pfamA_proteins_no_pfam_match.csv"
SUMMARY_REPORT = "pfamA_full_summary.txt"

CPUS = 4

# =========================
# CHECK INPUT FILES
# =========================
print("Step 1: Checking input files...")

for f in [HMM_FILE, FASTA_FILE]:
    if not Path(f).exists():
        raise FileNotFoundError(f"Missing input file: {f}")

# =========================
# CHECK HMMER INSTALLATION
# =========================
print("Step 2: Checking HMMER installation...")

try:
    subprocess.run(["hmmpress", "-h"], capture_output=True, text=True, check=True)
except FileNotFoundError:
    raise RuntimeError("HMMER not found! Run: conda install -c bioconda hmmer")

# =========================
# INDEX HMM DATABASE IF NEEDED
# =========================
print("Step 3: Checking HMM database index...")

if not Path(f"{HMM_FILE}.h3m").exists():
    print(f"Indexing {HMM_FILE} (this only happens once)...")
    subprocess.run(["hmmpress", HMM_FILE], check=True)

# =========================
# RUN HMMSCAN: TABLE OUTPUT
# =========================
print("Step 4: Running hmmscan (table pass)...")

subprocess.run([
    "hmmscan",
    "--cpu", str(CPUS),
    "--cut_ga",
    "--tblout", TBL_OUT,
    HMM_FILE,
    FASTA_FILE
], check=True)

# =========================
# RUN HMMSCAN: ALIGNMENT OUTPUT
# =========================
print("Step 5: Running hmmscan (alignment pass)...")

with open(ALIGNMENT_OUT, "w") as out:
    subprocess.run([
        "hmmscan",
        "--cpu", str(CPUS),
        "--cut_ga",
        HMM_FILE,
        FASTA_FILE
    ], stdout=out, check=True)

# =========================
# READ ALL PROTEIN IDS FROM FASTA
# =========================
print("Step 6: Reading protein IDs from FASTA...")

all_proteins_set = set()
with open(FASTA_FILE, "r") as f:
    for line in f:
        if line.startswith(">"):
            protein_id = line[1:].strip().split()[0]
            all_proteins_set.add(protein_id)

# =========================
# PARSE TBL OUT
# =========================
print("Step 7: Parsing hmmscan table output...")

rows = []
with open(TBL_OUT, "r") as f:
    for line in f:
        if line.startswith("#") or not line.strip():
            continue

        parts = line.split()
        if len(parts) < 18:
            continue

        rows.append({
            "protein": parts[2],
            "pfam_name": parts[0],
            "pfam_accession": parts[1].split(".")[0],
            "evalue": float(parts[4]),
            "score": float(parts[5])
        })

df = pd.DataFrame(rows)

# =========================
# NO MATCH CASE
# =========================
if df.empty:
    print("No Pfam matches found.")

    pd.DataFrame(columns=["protein", "pfam_name", "pfam_accession", "evalue", "score"]).to_csv(PFAM_FULL_RESULTS, index=False)
    pd.DataFrame(columns=["protein", "num_pfams", "pfam_accession", "pfam_name"]).to_csv(PROTEINS_SINGLE, index=False)
    pd.DataFrame(columns=["protein", "num_pfams", "pfam_accession", "pfam_name"]).to_csv(PROTEINS_MULTI, index=False)
    pd.DataFrame(columns=["pfam_accession", "pfam_name", "num_proteins"]).to_csv(PFAMS_MULTI_PROTEINS, index=False)
    pd.DataFrame(columns=["pfam_accession", "pfam_name", "num_proteins"]).to_csv(PFAMS_SINGLE_PROTEIN, index=False)
    pd.DataFrame({"protein": sorted(all_proteins_set)}).to_csv(PROTEINS_NO_MATCH, index=False)

    summary_lines = [
        "Pfam-A Full Scan Summary",
        "========================",
        f"Total proteins: {len(all_proteins_set)}",
        "Proteins with at least one Pfam match: 0",
        f"Proteins without Pfam match: {len(all_proteins_set)}",
        "Total protein-Pfam associations: 0",
        "Unique Pfam domains detected: 0",
        "Proteins with exactly one Pfam: 0",
        "Proteins with multiple Pfams: 0",
        "Pfams found in a single protein: 0",
        "Pfams found in multiple proteins: 0",
        "",
        "Files generated:",
        f"- {PFAM_FULL_RESULTS}",
        f"- {PROTEINS_SINGLE}",
        f"- {PROTEINS_MULTI}",
        f"- {PFAMS_MULTI_PROTEINS}",
        f"- {PFAMS_SINGLE_PROTEIN}",
        f"- {PROTEINS_NO_MATCH}",
        f"- {SUMMARY_REPORT}",
        f"- {TBL_OUT}",
        f"- {ALIGNMENT_OUT}",
    ]

    with open(SUMMARY_REPORT, "w") as f:
        f.write("\n".join(summary_lines))

    print(f"Summary written to {SUMMARY_REPORT}")
    print("\n".join(summary_lines))

# =========================
# MATCH CASE
# =========================
else:
    print("Step 8: Processing results into CSV files...")

    # Remove duplicate protein-Pfam pairs
    pair_df = df.drop_duplicates(subset=["protein", "pfam_accession"]).copy()
    pair_df.to_csv(PFAM_FULL_RESULTS, index=False)

    # -------------------------
    # Proteins with single / multiple Pfams
    # -------------------------
    protein_counts = (
        pair_df.groupby("protein")["pfam_accession"]
        .nunique()
        .reset_index(name="num_pfams")
    )

    proteins_single = protein_counts[protein_counts["num_pfams"] == 1].merge(
        pair_df[["protein", "pfam_accession", "pfam_name"]],
        on="protein",
        how="left"
    ).drop_duplicates()

    proteins_multi = protein_counts[protein_counts["num_pfams"] > 1].merge(
        pair_df[["protein", "pfam_accession", "pfam_name"]],
        on="protein",
        how="left"
    ).drop_duplicates()

    proteins_single = proteins_single.sort_values(["protein", "pfam_accession"])
    proteins_multi = proteins_multi.sort_values(["protein", "pfam_accession"])

    proteins_single.to_csv(PROTEINS_SINGLE, index=False)
    proteins_multi.to_csv(PROTEINS_MULTI, index=False)

    # -------------------------
    # Pfams found in single / multiple proteins
    # -------------------------
    pfam_counts = (
        pair_df.groupby(["pfam_accession", "pfam_name"])["protein"]
        .nunique()
        .reset_index(name="num_proteins")
    )

    pfams_multi = pfam_counts[pfam_counts["num_proteins"] > 1].copy()
    pfams_single = pfam_counts[pfam_counts["num_proteins"] == 1].copy()

    pfams_multi = pfams_multi.sort_values(["num_proteins", "pfam_accession"], ascending=[False, True])
    pfams_single = pfams_single.sort_values(["pfam_accession"])

    pfams_multi.to_csv(PFAMS_MULTI_PROTEINS, index=False)
    pfams_single.to_csv(PFAMS_SINGLE_PROTEIN, index=False)

    # -------------------------
    # Proteins with no Pfam match
    # -------------------------
    matched_proteins = set(pair_df["protein"].unique())
    unmatched_proteins = sorted(all_proteins_set - matched_proteins)

    pd.DataFrame({"protein": unmatched_proteins}).to_csv(PROTEINS_NO_MATCH, index=False)

    # -------------------------
    # Summary statistics
    # -------------------------
    total_proteins = len(all_proteins_set)
    proteins_with_match = len(matched_proteins)
    proteins_without_match = len(unmatched_proteins)
    total_associations = len(pair_df)
    unique_pfams = pair_df["pfam_accession"].nunique()
    proteins_exactly_one = (protein_counts["num_pfams"] == 1).sum()
    proteins_multiple = (protein_counts["num_pfams"] > 1).sum()
    pfams_single_count = (pfam_counts["num_proteins"] == 1).sum()
    pfams_multiple_count = (pfam_counts["num_proteins"] > 1).sum()

    summary_lines = [
        "Pfam-A Full Scan Summary",
        "========================",
        f"Total proteins: {total_proteins}",
        f"Proteins with at least one Pfam match: {proteins_with_match}",
        f"Proteins without Pfam match: {proteins_without_match}",
        f"Total protein-Pfam associations: {total_associations}",
        f"Unique Pfam domains detected: {unique_pfams}",
        f"Proteins with exactly one Pfam: {proteins_exactly_one}",
        f"Proteins with multiple Pfams: {proteins_multiple}",
        f"Pfams found in a single protein: {pfams_single_count}",
        f"Pfams found in multiple proteins: {pfams_multiple_count}",
        "",
        "Files generated:",
        f"- {PFAM_FULL_RESULTS}",
        f"- {PROTEINS_SINGLE}",
        f"- {PROTEINS_MULTI}",
        f"- {PFAMS_MULTI_PROTEINS}",
        f"- {PFAMS_SINGLE_PROTEIN}",
        f"- {PROTEINS_NO_MATCH}",
        f"- {SUMMARY_REPORT}",
        f"- {TBL_OUT}",
        f"- {ALIGNMENT_OUT}",
    ]

    with open(SUMMARY_REPORT, "w") as f:
        f.write("\n".join(summary_lines))

    print(f"Success! Summary written to {SUMMARY_REPORT}")
    print()
    print("\n".join(summary_lines))