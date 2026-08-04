import pandas as pd

# =========================
# INPUT
# =========================
INPUT_FILE = "pfamA_full_results.csv"

# =========================
# READ RESULTS
# =========================
df = pd.read_csv(INPUT_FILE)

required_cols = {"protein", "pfam_accession", "pfam_name"}
missing = required_cols - set(df.columns)
if missing:
    raise ValueError(f"Missing required columns in {INPUT_FILE}: {missing}")

# =========================
# FILTER DUF MATCHES
# =========================
duf_df = df[df["pfam_name"].astype(str).str.contains(r"DUF", case=False, na=False)].copy()

# Optional: remove duplicate protein-domain pairs
duf_df = duf_df.drop_duplicates(subset=["protein", "pfam_accession"])

# =========================
# SAVE FULL DUF MATCH TABLE
# =========================
duf_df.to_csv("e.coli_k-12_pfama_all_DUF_matches.csv", index=False)

# =========================
# DUF COUNTS PER DOMAIN
# =========================
duf_domain_counts = (
    duf_df.groupby(["pfam_accession", "pfam_name"])["protein"]
    .nunique()
    .reset_index(name="num_proteins")
    .sort_values(["num_proteins", "pfam_accession"], ascending=[False, True])
)

duf_domain_counts.to_csv("e.coli_k-12_pfama_DUF_domain_counts.csv", index=False)

# =========================
# DUF COUNTS PER PROTEIN
# =========================
protein_duf_counts = (
    duf_df.groupby("protein")["pfam_accession"]
    .nunique()
    .reset_index(name="num_DUFs")
    .sort_values(["num_DUFs", "protein"], ascending=[False, True])
)

protein_duf_counts.to_csv("e.coli_k-12_pfama_proteins_with_DUFs.csv", index=False)

# =========================
# SUMMARY NUMBERS
# =========================
total_matches = len(df)
total_unique_pfams = df["pfam_accession"].nunique()
total_proteins_with_any_pfam = df["protein"].nunique()

total_duf_matches = len(duf_df)
total_unique_dufs = duf_df["pfam_accession"].nunique()
total_duf_names = duf_df["pfam_name"].nunique()
proteins_with_dufs = duf_df["protein"].nunique()

# proteins with exactly one DUF vs multiple DUFs
proteins_one_duf = (protein_duf_counts["num_DUFs"] == 1).sum()
proteins_multiple_dufs = (protein_duf_counts["num_DUFs"] > 1).sum()

# DUFs found in exactly one protein vs multiple proteins
dufs_one_protein = (duf_domain_counts["num_proteins"] == 1).sum()
dufs_multiple_proteins = (duf_domain_counts["num_proteins"] > 1).sum()

# =========================
# TOP TABLES FOR SUMMARY
# =========================
top_dufs_txt = (
    duf_domain_counts.head(10).to_string(index=False)
    if not duf_domain_counts.empty else "No DUF matches found"
)

top_proteins_txt = (
    protein_duf_counts.head(10).to_string(index=False)
    if not protein_duf_counts.empty else "No proteins with DUFs found"
)

# =========================
# SUMMARY REPORT
# =========================
summary = f"""
e.coli_k-12 FULL PFAM-A DUF MATCH ANALYSIS
===================================

Input file:
- {INPUT_FILE}

Overall Pfam-A annotation:
- Total Pfam matches: {total_matches}
- Unique Pfam domains detected: {total_unique_pfams}
- Proteins with at least one Pfam match: {total_proteins_with_any_pfam}

DUF-specific results:
- Total DUF matches: {total_duf_matches}
- Unique DUF domains detected: {total_unique_dufs}
- Unique DUF names detected: {total_duf_names}
- Proteins containing at least one DUF: {proteins_with_dufs}

Protein-level DUF distribution:
- Proteins with exactly one DUF: {proteins_one_duf}
- Proteins with multiple DUFs: {proteins_multiple_dufs}

DUF-level protein distribution:
- DUFs found in exactly one protein: {dufs_one_protein}
- DUFs found in multiple proteins: {dufs_multiple_proteins}

Top DUF domains by number of proteins:
{top_dufs_txt}

Top proteins by number of DUF domains:
{top_proteins_txt}

Files generated:
- e.coli_k-12_pfama_all_DUF_matches.csv
- e.coli_k-12_pfama_DUF_domain_counts.csv
- e.coli_k-12_pfama_proteins_with_DUFs.csv
- e.coli_k-12_pfama_DUF_summary.txt
"""

with open("e.coli_k-12_pfama_DUF_summary.txt", "w") as f:
    f.write(summary)

print(summary)
