library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(scales)

# =========================
# INPUT FILES
# =========================
ESSENTIAL_FILE <- "essential102_pfam_full_results.csv"
FASTAAI_FILE   <- "fastaai122_pfam_full_results.csv"
PFAM500_FILE   <- "pfam_500_full_results.csv"
PFAMFULL_FILE  <- "pfamA_full_results.csv"
HYP_FILE       <- "all_hypothetical_proteins.csv"
DUF_FILE       <- "e.coli_k-12_pfama_all_DUF_matches.csv"

TOTAL_PROTEINS <- 4300

# =========================
# LOAD DATA
# =========================
essential <- read_csv(ESSENTIAL_FILE, show_col_types = FALSE)
fastaai   <- read_csv(FASTAAI_FILE, show_col_types = FALSE)
pfam500   <- read_csv(PFAM500_FILE, show_col_types = FALSE)
pfam_full <- read_csv(PFAMFULL_FILE, show_col_types = FALSE)
hyp       <- read_csv(HYP_FILE, show_col_types = FALSE)
duf       <- read_csv(DUF_FILE, show_col_types = FALSE)

# =========================
# HELPER FUNCTION
# =========================
make_summary <- function(df, total_proteins, set_name) {
  pair_df <- df %>%
    distinct(protein, pfam_accession, .keep_all = TRUE)

  protein_counts <- pair_df %>%
    count(protein, name = "num_pfams")

  data.frame(
    set = set_name,
    total_proteins = total_proteins,
    matched_proteins = n_distinct(pair_df$protein),
    unmatched_proteins = total_proteins - n_distinct(pair_df$protein),
    unique_pfams = n_distinct(pair_df$pfam_accession),
    single_domain_proteins = sum(protein_counts$num_pfams == 1),
    multi_domain_proteins  = sum(protein_counts$num_pfams > 1)
  )
}

summary_df <- bind_rows(
  make_summary(essential, TOTAL_PROTEINS, "Essential 102"),
  make_summary(fastaai,   TOTAL_PROTEINS, "FastAAI 122"),
  make_summary(pfam500,   TOTAL_PROTEINS, "Pfam 500"),
  make_summary(pfam_full, TOTAL_PROTEINS, "Pfam-A full")
)

summary_df$set <- factor(
  summary_df$set,
  levels = c("Essential 102", "FastAAI 122", "Pfam 500", "Pfam-A full")
)

# =========================
# COMMON THEME
# =========================
theme_clean <- theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 15),
    plot.subtitle = element_text(size = 10),
    axis.text.x = element_text(angle = 25, hjust = 1, size = 10),
    axis.text.y = element_text(size = 10),
    axis.title = element_text(size = 12),
    legend.position = "top",
    legend.title = element_blank(),
    legend.text = element_text(size = 10),
    panel.grid.major = element_line(color = "grey85"),
    panel.grid.minor = element_line(color = "grey92")
  )

# =========================
# PANEL A: COVERAGE (ROBUST LABEL FIX)
# =========================
coverage_df <- summary_df %>%
  select(set, matched_proteins, unmatched_proteins) %>%
  pivot_longer(
    cols = c(matched_proteins, unmatched_proteins),
    names_to = "category",
    values_to = "count"
  ) %>%
  group_by(set) %>%
  mutate(
    total = sum(count),
    frac = count / total,
    label = paste0(count, "\n(", round(frac * 100), "%)")
  ) %>%
  ungroup()

coverage_labels <- coverage_df %>%
  select(set, category, count, frac, label) %>%
  pivot_wider(names_from = category, values_from = c(count, frac, label)) %>%
  mutate(
    y_unmatched = frac_unmatched_proteins / 2,
    y_matched_inside = frac_unmatched_proteins + frac_matched_proteins / 2,
    y_matched_outside = 1.02
  )

small_cutoff_A <- 0.06

pA <- ggplot(coverage_df, aes(x = set, y = frac, fill = category)) +
  geom_col(width = 0.72, color = "black") +

  # lower segment label always inside
  geom_text(
    data = coverage_labels,
    aes(x = set, y = y_unmatched, label = label_unmatched_proteins),
    inherit.aes = FALSE,
    size = 3.7,
    fontface = "bold",
    lineheight = 0.9
  ) +

  # top segment labels INSIDE if large enough
  geom_text(
    data = coverage_labels %>% filter(frac_matched_proteins >= small_cutoff_A),
    aes(x = set, y = y_matched_inside, label = label_matched_proteins),
    inherit.aes = FALSE,
    size = 3.7,
    fontface = "bold",
    lineheight = 0.9
  ) +

  # top segment labels OUTSIDE if too small
  geom_text(
    data = coverage_labels %>% filter(frac_matched_proteins < small_cutoff_A),
    aes(x = set, y = y_matched_outside, label = label_matched_proteins),
    inherit.aes = FALSE,
    size = 3.7,
    fontface = "bold",
    lineheight = 0.9,
    vjust = 0
  ) +

  scale_y_continuous(
    labels = scales::percent_format(),
    limits = c(0, 1.10),
    expand = c(0, 0)
  ) +
  coord_cartesian(clip = "off") +
  scale_fill_manual(
    values = c(
      "matched_proteins" = "#5B8CC0",
      "unmatched_proteins" = "#D97941"
    ),
    labels = c(
      "matched_proteins" = "With Pfam",
      "unmatched_proteins" = "Without Pfam"
    )
  ) +
  labs(
    title = "A. Protein coverage across Pfam sets",
    subtitle = "Fraction of E. coli K-12 proteins matched by each Pfam set",
    x = NULL,
    y = "Fraction of proteins",
    fill = NULL
  ) +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 15),
    plot.subtitle = element_text(size = 10),
    axis.text.x = element_text(angle = 25, hjust = 1, size = 10),
    legend.position = "top",
    plot.margin = margin(10, 10, 10, 10)
  )
# =========================
# PANEL B
# =========================
pB <- ggplot(summary_df, aes(x = set, y = unique_pfams)) +
  geom_col(
    width = 0.68,
    fill = c("#4F79A7", "#59A96A", "#C44E52", "#8172B2"),
    color = "black"
  ) +
  geom_text(
    aes(label = unique_pfams),
    vjust = -0.35,
    size = 3.8,
    fontface = "bold"
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.06))) +
  labs(
    title = "B. Unique Pfam domains detected",
    subtitle = "Functional diversity captured by each Pfam set",
    x = NULL,
    y = "Number of unique Pfam domains"
  ) +
  theme_clean

# =========================
# PANEL C: ARCHITECTURE (ROBUST LABEL FIX)
# =========================
arch_df <- summary_df %>%
  select(set, single_domain_proteins, multi_domain_proteins) %>%
  pivot_longer(
    cols = c(single_domain_proteins, multi_domain_proteins),
    names_to = "category",
    values_to = "count"
  )

arch_labels <- arch_df %>%
  select(set, category, count) %>%
  pivot_wider(names_from = category, values_from = count) %>%
  mutate(
    total = single_domain_proteins + multi_domain_proteins,
    y_single = single_domain_proteins / 2,
    y_multi_inside = single_domain_proteins + multi_domain_proteins / 2,
    y_multi_outside = total + max(total) * 0.03
  )

small_cutoff_C <- 120

pC <- ggplot(arch_df, aes(x = set, y = count, fill = category)) +
  geom_col(width = 0.72, color = "black") +

  # bottom segment labels
  geom_text(
    data = arch_labels,
    aes(x = set, y = y_single, label = single_domain_proteins),
    inherit.aes = FALSE,
    size = 3.7,
    fontface = "bold"
  ) +

  # top segment labels inside if large
  geom_text(
    data = arch_labels %>% filter(multi_domain_proteins >= small_cutoff_C),
    aes(x = set, y = y_multi_inside, label = multi_domain_proteins),
    inherit.aes = FALSE,
    size = 3.7,
    fontface = "bold"
  ) +

  # top segment labels outside if tiny
  geom_text(
    data = arch_labels %>% filter(multi_domain_proteins < small_cutoff_C),
    aes(x = set, y = y_multi_outside, label = multi_domain_proteins),
    inherit.aes = FALSE,
    size = 3.7,
    fontface = "bold",
    vjust = 0
  ) +

  scale_fill_manual(
    values = c(
      "multi_domain_proteins" = "#E5B300",
      "single_domain_proteins" = "#8DA0CB"
    ),
    labels = c(
      "multi_domain_proteins" = "Multi-domain proteins",
      "single_domain_proteins" = "Single-domain proteins"
    )
  ) +
  scale_y_continuous(
    limits = c(0, max(arch_labels$total) * 1.10),
    expand = c(0, 0)
  ) +
  coord_cartesian(clip = "off") +
  labs(
    title = "C. Protein domain architecture across Pfam sets",
    subtitle = "Single-domain versus multi-domain proteins",
    x = NULL,
    y = "Number of proteins",
    fill = NULL
  ) +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 15),
    plot.subtitle = element_text(size = 10),
    axis.text.x = element_text(angle = 25, hjust = 1, size = 10),
    legend.position = "top",
    plot.margin = margin(10, 10, 10, 10)
  )
# =========================
# PANEL D
# =========================
pfam_full_pairs <- pfam_full %>%
  distinct(protein, pfam_accession)

pfam_annotated <- n_distinct(pfam_full_pairs$protein)
no_pfam <- TOTAL_PROTEINS - pfam_annotated
duf_proteins <- n_distinct(duf$protein)

landscape_df <- data.frame(
  category = factor(
    c("Total proteins", "Pfam-A annotated", "DUF proteins", "No Pfam match"),
    levels = c("Total proteins", "Pfam-A annotated", "DUF proteins", "No Pfam match")
  ),
  count = c(TOTAL_PROTEINS, pfam_annotated, duf_proteins, no_pfam)
)

pD <- ggplot(landscape_df, aes(x = category, y = count)) +
  geom_col(
    width = 0.68,
    fill = c("#5B8CC0", "#A8A8A8", "#8FB0D9", "#D97941"),
    color = "black"
  ) +
  geom_text(
    aes(label = count),
    vjust = -0.35,
    size = 3.8,
    fontface = "bold"
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(
    title = "D. Protein annotation landscape of E. coli K-12",
    subtitle = "Global Pfam-A annotation and DUF overview",
    x = NULL,
    y = "Number of proteins"
  ) +
  theme_clean

# =========================
# PANEL E
# =========================
hyp_with_pfam <- pfam_full_pairs %>%
  filter(protein %in% hyp$protein) %>%
  distinct(protein) %>%
  nrow()

hyp_with_duf <- duf %>%
  filter(protein %in% hyp$protein) %>%
  distinct(protein) %>%
  nrow()

hyp_total <- nrow(hyp)
hyp_no_pfam <- hyp_total - hyp_with_pfam

hyp_df <- data.frame(
  category = factor(
    c("Total hypothetical", "Pfam-A annotated", "DUF proteins", "No Pfam match"),
    levels = c("Total hypothetical", "Pfam-A annotated", "DUF proteins", "No Pfam match")
  ),
  count = c(hyp_total, hyp_with_pfam, hyp_with_duf, hyp_no_pfam)
)

pE <- ggplot(hyp_df, aes(x = category, y = count)) +
  geom_col(
    width = 0.68,
    fill = c("#5B8CC0", "#A8A8A8", "#8FB0D9", "#D97941"),
    color = "black"
  ) +
  geom_text(
    aes(label = count),
    vjust = -0.35,
    size = 3.8,
    fontface = "bold"
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.06))) +
  labs(
    title = paste0("E. Hypothetical protein annotation (n = ", hyp_total, ")"),
    subtitle = "Pfam-A and DUF characterization of hypothetical proteins",
    x = NULL,
    y = "Number of proteins"
  ) +
  theme_clean +
  theme(
    axis.text.x = element_text(angle = 20, hjust = 1, size = 10)
  )

# =========================
# COMBINE
# =========================
final_plot <- (pA + pB) / (pC + pD) / pE +
  plot_annotation(
    title = "Protein domain annotation in E. coli K-12",
    subtitle = "Comparison of targeted Pfam sets, global Pfam-A annotation, DUF content, and hypothetical protein characterization",
    theme = theme(
      plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 11, hjust = 0.5)
    )
  ) &
  theme(
    plot.margin = margin(8, 10, 8, 10)
  )

print(final_plot)

ggsave("EcoliK12_Pfam_collage_fixed.png", final_plot, width = 15, height = 18, dpi = 300)
ggsave("EcoliK12_Pfam_collage_fixed.pdf", final_plot, width = 15, height = 18)
