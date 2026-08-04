library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)

# =========================
# FILES
# =========================
pfam_full <- read_csv("pfamA_full_results.csv", show_col_types = FALSE)
essential <- read_csv("essential102_pfam_full_results.csv", show_col_types = FALSE)
fastaai   <- read_csv("fastaai122_pfam_full_results.csv", show_col_types = FALSE)
pfam500   <- read_csv("pfam_500_full_results.csv", show_col_types = FALSE)
hyp       <- read_csv("all_hypothetical_proteins.csv", show_col_types = FALSE)
duf       <- read_csv("e.coli_k-12_pfama_all_DUF_matches.csv", show_col_types = FALSE)

TOTAL_PROTEINS <- 4300   # change if needed

# =========================
# FUNCTION
# =========================
get_summary <- function(df, name) {
  pairs <- df %>% distinct(protein, pfam_accession)

  data.frame(
    Category = name,
    Proteins_matched = n_distinct(pairs$protein),
    Detected_Pfams   = n_distinct(pairs$pfam_accession),
    DUF_proteins     = NA,
    No_Pfam_match    = TOTAL_PROTEINS - n_distinct(pairs$protein)
  )
}

# =========================
# BUILD SUMMARY
# =========================
summary_df <- bind_rows(
  get_summary(essential, "Essential 102"),
  get_summary(fastaai,   "FastAAI 122"),
  get_summary(pfam500,   "Pfam 500")
)

# Global Pfam-A
pfam_pairs <- pfam_full %>% distinct(protein, pfam_accession)

global <- data.frame(
  Category = "Global Pfam-A",
  Proteins_matched = n_distinct(pfam_pairs$protein),
  Detected_Pfams   = n_distinct(pfam_pairs$pfam_accession),
  DUF_proteins     = n_distinct(duf$protein),
  No_Pfam_match    = TOTAL_PROTEINS - n_distinct(pfam_pairs$protein)
)

# Hypothetical
hyp_with_pfam <- pfam_pairs %>%
  filter(protein %in% hyp$protein) %>%
  distinct(protein) %>%
  nrow()

hyp_with_duf <- duf %>%
  filter(protein %in% hyp$protein) %>%
  distinct(protein) %>%
  nrow()

hyp_sum <- data.frame(
  Category = "Hypothetical",
  Proteins_matched = hyp_with_pfam,
  Detected_Pfams   = pfam_full %>%
    filter(protein %in% hyp$protein) %>%
    distinct(pfam_accession) %>%
    nrow(),
  DUF_proteins     = hyp_with_duf,
  No_Pfam_match    = nrow(hyp) - hyp_with_pfam
)

summary_df <- bind_rows(summary_df, global, hyp_sum)

# =========================
# FORMAT FOR PLOTTING
# =========================
plot_df <- summary_df %>%
  pivot_longer(
    cols = c(Proteins_matched, Detected_Pfams, DUF_proteins, No_Pfam_match),
    names_to = "Metric",
    values_to = "Count"
  ) %>%
  filter(!is.na(Count))

plot_df$Category <- factor(
  plot_df$Category,
  levels = c("Essential 102", "FastAAI 122", "Pfam 500", "Global Pfam-A", "Hypothetical")
)

plot_df$Metric <- factor(
  plot_df$Metric,
  levels = c("Proteins_matched", "Detected_Pfams", "DUF_proteins", "No_Pfam_match"),
  labels = c("Proteins matched", "Detected Pfams", "DUF proteins", "No Pfam match")
)

# =========================
# PLOT (IDENTICAL STYLE)
# =========================
p <- ggplot(plot_df, aes(x = Category, y = Count, fill = Metric)) +
  geom_col(
    position = position_dodge(width = 0.75),
    width = 0.65,
    color = "black",
    linewidth = 0.6
  ) +
  geom_text(
    aes(label = Count),
    position = position_dodge(width = 0.75),
    vjust = -0.3,
    size = 5,
    fontface = "bold"
  ) +
  scale_fill_manual(
    values = c(
      "Proteins matched" = "#4F79A7",
      "Detected Pfams"   = "#59A14F",
      "DUF proteins"     = "#B07AA1",
      "No Pfam match"    = "#E15759"
    )
  ) +
  labs(
    title = "E. coli K-12 protein domain annotation summary",
    x = NULL,
    y = "Count",
    fill = NULL
  ) +
  theme_minimal(base_size = 18) +
  theme(
    plot.title = element_text(face = "bold", size = 22, hjust = 0.5),
    axis.text.x = element_text(angle = 15, hjust = 1, face = "bold"),
    axis.title.y = element_text(face = "bold"),
    legend.position = "top",
    panel.grid.major = element_line(color = "grey80"),
    panel.grid.minor = element_line(color = "grey90")
  )

print(p)

ggsave("Ecoli_K12_final_plot.png", p, width = 16, height = 9, dpi = 300)
ggsave("Ecoli_K12_final_plot.pdf", p, width = 16, height = 9)
