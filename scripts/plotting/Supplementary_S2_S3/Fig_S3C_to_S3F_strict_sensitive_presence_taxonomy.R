#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(cowplot)
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(stringr)
  library(tidyr)
})

repo_root <- Sys.getenv("PAPER_SUBMISSION_CODE", unset = normalizePath(getwd(), winslash = "/", mustWork = TRUE))
project_dir <- file.path(repo_root, "data", "external", "supplementary_s2_s3", "human_hq_mag_cf_matrix")
analysis_dir <- file.path(project_dir, "main_strict_sensitive_result", "cf_presence_absence_analysis")
tables_dir <- file.path(analysis_dir, "tables")
fig_dir <- file.path(repo_root, "outputs", "Supplementary_S2_S3")
plot_table_dir <- file.path(fig_dir, "plot_tables")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_table_dir, recursive = TRUE, showWarnings = FALSE)

source(file.path(repo_root, "scripts", "plotting", "Figure3", "Figure3_color_config.R"))

genus_order <- c("Staphylococcus", "Corynebacterium", "Cutibacterium")
min_species_mags <- 20

summary_df <- read_tsv(file.path(tables_dir, "cf_presence_taxonomy_summary_for_plot.tsv"), show_col_types = FALSE) %>%
  mutate(
    eCM_label = str_replace(eCM, "module", "eCM"),
    CF_label = str_replace(CF_id, "^CF", "eCF"),
    module_fill = unname(full_module_colors[eCM]),
    present_for_plot = prevalence > 0,
    copy_bin = case_when(
      !present_for_plot ~ NA_character_,
      mean_copy_present <= 1 ~ "1",
      mean_copy_present <= 1.25 ~ ">1-1.25",
      mean_copy_present <= 1.5 ~ ">1.25-1.5",
      mean_copy_present <= 2 ~ ">1.5-2",
      mean_copy_present <= 3 ~ ">2-3",
      TRUE ~ ">3"
    ),
    copy_bin = factor(copy_bin, levels = c("1", ">1-1.25", ">1.25-1.5", ">1.5-2", ">2-3", ">3")),
    significant_global = !is.na(global_FDR) & global_FDR < 0.05
  )

cf_order <- summary_df %>%
  distinct(CF_id, CF_label, eCM, eCM_label, eCM_num, module_fill) %>%
  mutate(
    cf_major = as.integer(str_extract(CF_id, "(?<=CF)\\d+")),
    cf_minor = as.integer(str_extract(CF_id, "(?<=_)\\d+")),
    cf_minor = if_else(is.na(cf_minor), -1L, cf_minor)
  ) %>%
  arrange(eCM_num, cf_major, cf_minor, CF_id) %>%
  mutate(x_pos = row_number())

copy_bin_colors <- c(
  "1" = "#D9D9D9",
  ">1-1.25" = "#FEE0D2",
  ">1.25-1.5" = "#FCBBA1",
  ">1.5-2" = "#FC9272",
  ">2-3" = "#FB6A4A",
  ">3" = "#A50F15"
)

module_bar <- cf_order %>%
  transmute(CF_id, x_pos, eCM_label, module_fill)

write_tsv(summary_df, file.path(plot_table_dir, "cf_presence_taxonomy_summary_for_plot_with_bins.tsv"))

make_sample_bar <- function(layout_df, y_limits, width_title = "Number of MAGs") {
  n_max <- max(layout_df$n_mags, na.rm = TRUE)
  ggplot(layout_df, aes(x = n_mags, y = y_pos)) +
    geom_col(width = 0.58, fill = "#4B5563", color = "#374151", linewidth = 0.15, orientation = "y") +
    geom_text(aes(label = n_mags), hjust = -0.12, size = 2.4, color = "black") +
    scale_y_continuous(limits = y_limits, breaks = layout_df$y_pos, labels = NULL, expand = c(0, 0)) +
    scale_x_continuous(limits = c(0, n_max * 1.35), breaks = scales::breaks_pretty(n = 3), expand = c(0, 0)) +
    labs(x = width_title, y = NULL, title = "Sample size") +
    theme_minimal(base_size = 8) +
    theme(
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.text.x = element_text(size = 6, color = "black"),
      axis.title.x = element_text(size = 7, color = "black", margin = margin(t = 2)),
      plot.title = element_text(size = 8, face = "bold", hjust = 0.5),
      plot.margin = margin(24, 2, 43, 0)
    )
}

plot_presence_bubble <- function(plot_df, layout_df, sig_df, separators, title, subtitle, file_base, width, height) {
  y_top <- max(layout_df$y_pos) + 1
  y_limits <- c(min(layout_df$y_pos) - 0.7, y_top + 0.9)
  annot_df <- module_bar %>% mutate(y_pos = y_top)
  p_main <- ggplot() +
    geom_hline(yintercept = separators, color = "grey70", linewidth = 0.25) +
    geom_tile(
      data = annot_df,
      aes(x = x_pos, y = y_pos, fill = module_fill),
      width = 0.96, height = 0.55, color = "white", linewidth = 0.25
    ) +
    geom_text(data = annot_df, aes(x = x_pos, y = y_pos, label = eCM_label), size = 1.85, color = "black") +
    geom_text(
      data = sig_df,
      aes(x = x_pos, y = y_top + 0.5, label = "*"),
      size = 3.1, color = "black"
    ) +
    geom_point(
      data = plot_df %>% filter(present_for_plot),
      aes(x = x_pos, y = y_pos, size = prevalence, color = copy_bin),
      alpha = 0.95
    ) +
    scale_fill_identity() +
    scale_color_manual(
      name = "Mean copy\npresent MAGs",
      values = copy_bin_colors,
      limits = names(copy_bin_colors),
      drop = FALSE
    ) +
    scale_size_continuous(
      name = "eCF-positive MAG prevalence",
      range = c(0.5, 7),
      limits = c(0, 1),
      breaks = c(0.25, 0.5, 0.75, 1),
      labels = scales::percent_format(accuracy = 1)
    ) +
    scale_x_continuous(breaks = cf_order$x_pos, labels = cf_order$CF_label, expand = expansion(mult = c(0.005, 0.005))) +
    scale_y_continuous(limits = y_limits, breaks = layout_df$y_pos, labels = layout_df$taxon_label, expand = c(0, 0)) +
    labs(x = NULL, y = NULL, title = title, subtitle = subtitle) +
    theme_minimal(base_size = 8) +
    theme(
      panel.grid.major = element_line(color = "grey91", linewidth = 0.25),
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 7, color = "black"),
      axis.text.y = element_text(size = 8, color = "black"),
      legend.position = "bottom",
      legend.box = "horizontal",
      plot.title = element_text(size = 11, face = "bold"),
      plot.subtitle = element_text(size = 8, color = "grey30"),
      plot.margin = margin(5, 5, 5, 5)
    )
  p_n <- make_sample_bar(layout_df, y_limits)
  p <- cowplot::plot_grid(p_main, p_n, nrow = 1, rel_widths = c(1, 0.18), align = "h", axis = "tb")
  ggsave(file.path(fig_dir, paste0(file_base, ".png")), p, width = width, height = height, dpi = 300)
  ggsave(file.path(fig_dir, paste0(file_base, ".pdf")), p, width = width, height = height)
}

genus_layout <- summary_df %>%
  filter(taxonomic_level == "genus") %>%
  distinct(taxon, n_mags) %>%
  mutate(
    taxon = factor(taxon, levels = genus_order),
    y_pos = length(genus_order) - match(as.character(taxon), genus_order) + 1,
    taxon_label = as.character(taxon)
  ) %>%
  arrange(y_pos)

genus_plot <- summary_df %>%
  filter(taxonomic_level == "genus") %>%
  left_join(genus_layout %>% select(taxon, y_pos), by = "taxon") %>%
  left_join(cf_order %>% select(CF_id, x_pos), by = "CF_id")

genus_sig <- genus_plot %>%
  filter(significant_global) %>%
  distinct(CF_id, x_pos)

plot_presence_bubble(
  genus_plot,
  genus_layout,
  genus_sig,
  numeric(0),
  "Strict-sensitive genus-level eCF presence and copy deployment",
  "Point size shows eCF-positive MAG prevalence; color shows mean copy number among eCF-positive MAGs. * marks genus-level FDR < 0.05.",
  "Fig_S3C_strict_sensitive_genus_eCF_presence_copy_deployment",
  13.8, 4.2
)

species_layout <- summary_df %>%
  filter(taxonomic_level == "species") %>%
  distinct(genus, taxon, n_mags) %>%
  mutate(genus = factor(genus, levels = genus_order)) %>%
  arrange(genus, desc(n_mags), taxon) %>%
  group_by(genus) %>%
  mutate(rank_in_genus = row_number()) %>%
  ungroup() %>%
  mutate(
    genus_index = as.integer(genus),
    y_raw = row_number() + (genus_index - 1) * 0.8,
    y_pos = max(y_raw) - y_raw + 1,
    taxon_label = taxon
  )

separator_df <- species_layout %>%
  arrange(y_pos) %>%
  mutate(next_genus = lead(genus), next_y = lead(y_pos)) %>%
  filter(!is.na(next_genus), next_genus != genus) %>%
  transmute(sep = (y_pos + next_y) / 2)

species_plot <- summary_df %>%
  filter(taxonomic_level == "species") %>%
  left_join(species_layout %>% select(taxon, y_pos), by = "taxon") %>%
  left_join(cf_order %>% select(CF_id, x_pos), by = "CF_id")

species_sig <- species_plot %>%
  filter(significant_global) %>%
  distinct(genus, CF_id, x_pos)

plot_presence_bubble(
  species_plot,
  species_layout,
  species_sig,
  separator_df$sep,
  paste0("Strict-sensitive within-genus species eCF deployment (n >= ", min_species_mags, " MAGs)"),
  "Species are grouped by genus with gaps; * marks within-genus global FDR < 0.05 for that genus-eCF test.",
  "Fig_S3D_strict_sensitive_within_genus_species_eCF_deployment",
  14.5, 8.2
)

genus_assoc <- read_tsv(file.path(tables_dir, "cf_presence_genus_global_association.tsv"), show_col_types = FALSE) %>%
  mutate(
    CF_label = str_replace(CF_id, "^CF", "eCF"),
    sig_class = case_when(FDR < 0.001 ~ "FDR < 0.001", FDR < 0.01 ~ "FDR < 0.01", FDR < 0.05 ~ "FDR < 0.05", TRUE ~ "FDR >= 0.05")
  )

species_assoc <- read_tsv(file.path(tables_dir, "cf_presence_species_within_genus_global_association.tsv"), show_col_types = FALSE) %>%
  mutate(
    CF_label = str_replace(CF_id, "^CF", "eCF"),
    sig_class = case_when(FDR < 0.001 ~ "FDR < 0.001", FDR < 0.01 ~ "FDR < 0.01", FDR < 0.05 ~ "FDR < 0.05", TRUE ~ "FDR >= 0.05")
  )

rank_genus <- genus_assoc %>%
  arrange(prevalence_range) %>%
  mutate(CF_label = factor(CF_label, levels = CF_label))

p_rank_genus <- ggplot(rank_genus, aes(x = prevalence_range, y = CF_label, color = sig_class)) +
  geom_segment(aes(x = 0, xend = prevalence_range, yend = CF_label), linewidth = 0.45) +
  geom_point(size = 2.2) +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1), expand = expansion(mult = c(0, 0.04))) +
  scale_color_manual(values = c("FDR < 0.001" = "#A50F15", "FDR < 0.01" = "#FB6A4A", "FDR < 0.05" = "#FC9272", "FDR >= 0.05" = "grey55"), drop = FALSE) +
  labs(x = "Max prevalence - min prevalence", y = NULL, color = NULL, title = "Genus-level eCF prevalence range") +
  theme_minimal(base_size = 8) +
  theme(panel.grid.major.y = element_blank(), axis.text.y = element_text(size = 7, color = "black"), plot.title = element_text(size = 11, face = "bold"))

rank_species <- species_assoc %>%
  group_by(genus) %>%
  slice_max(prevalence_range, n = 12, with_ties = FALSE) %>%
  ungroup() %>%
  arrange(genus, prevalence_range) %>%
  mutate(CF_label_rank = factor(paste(genus, CF_label, sep = " | "), levels = paste(genus, CF_label, sep = " | ")))

p_rank_species <- ggplot(rank_species, aes(x = prevalence_range, y = CF_label_rank, color = sig_class)) +
  geom_segment(aes(x = 0, xend = prevalence_range, yend = CF_label_rank), linewidth = 0.45) +
  geom_point(size = 2.2) +
  facet_grid(genus ~ ., scales = "free_y", space = "free_y") +
  scale_y_discrete(labels = function(x) sub("^.* \\| ", "", x)) +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1), expand = expansion(mult = c(0, 0.04))) +
  scale_color_manual(values = c("FDR < 0.001" = "#A50F15", "FDR < 0.01" = "#FB6A4A", "FDR < 0.05" = "#FC9272", "FDR >= 0.05" = "grey55"), drop = FALSE) +
  labs(x = "Within-genus max species prevalence - min species prevalence", y = NULL, color = NULL, title = "Species-level eCF prevalence range within each genus") +
  theme_minimal(base_size = 8) +
  theme(panel.grid.major.y = element_blank(), axis.text.y = element_text(size = 7, color = "black"), strip.text.y = element_text(face = "bold"), plot.title = element_text(size = 11, face = "bold"))

p_range <- cowplot::plot_grid(p_rank_genus, p_rank_species, nrow = 1, rel_widths = c(0.9, 1.15))
ggsave(file.path(fig_dir, "Fig_S3E_strict_sensitive_eCF_prevalence_range_rank.png"), p_range, width = 13, height = 7.2, dpi = 300)
ggsave(file.path(fig_dir, "Fig_S3E_strict_sensitive_eCF_prevalence_range_rank.pdf"), p_range, width = 13, height = 7.2)

pairwise <- read_tsv(file.path(tables_dir, "cf_presence_species_within_genus_pairwise_association.tsv"), show_col_types = FALSE)
if (nrow(pairwise) > 0) {
  pair_plot <- pairwise %>%
    mutate(
      CF_label = str_replace(CF_id, "^CF", "eCF"),
      pair_label = paste(species_a, species_b, sep = " vs "),
      pair_label = str_replace_all(pair_label, paste0(genus, " "), ""),
      sig_mark = if_else(FDR < 0.05, "*", "")
    )
  pair_levels <- pair_plot %>%
    group_by(genus, pair_label) %>%
    summarise(max_abs = max(abs_prevalence_diff, na.rm = TRUE), .groups = "drop") %>%
    arrange(genus, desc(max_abs), pair_label) %>%
    pull(pair_label) %>%
    unique()
  cf_levels <- pair_plot %>%
    group_by(CF_label) %>%
    summarise(max_abs = max(abs_prevalence_diff, na.rm = TRUE), .groups = "drop") %>%
    arrange(max_abs) %>%
    pull(CF_label)
  pair_plot <- pair_plot %>%
    mutate(pair_label = factor(pair_label, levels = pair_levels), CF_label = factor(CF_label, levels = cf_levels))
  write_tsv(pair_plot, file.path(plot_table_dir, "cf_presence_species_pairwise_heatmap_plot_data.tsv"))
  p_pair <- ggplot(pair_plot, aes(x = pair_label, y = CF_label, fill = prevalence_diff_a_minus_b)) +
    geom_tile(color = "white", linewidth = 0.18) +
    geom_text(aes(label = sig_mark), size = 2.5, color = "black") +
    facet_grid(. ~ genus, scales = "free_x", space = "free_x") +
    scale_fill_gradient2(
      name = "Prevalence diff\nspecies A - B",
      low = "#2166AC",
      mid = "#F7F7F7",
      high = "#B2182B",
      midpoint = 0,
      labels = scales::percent_format(accuracy = 1)
    ) +
    labs(x = NULL, y = NULL, title = "Pairwise species differences for within-genus significant eCFs") +
    theme_minimal(base_size = 8) +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 6.2, color = "black"),
      axis.text.y = element_text(size = 7, color = "black"),
      strip.text.x = element_text(face = "bold"),
      legend.position = "bottom",
      plot.title = element_text(size = 11, face = "bold")
    )
  ggsave(file.path(fig_dir, "Fig_S3F_strict_sensitive_within_genus_pairwise_species_heatmap.png"), p_pair, width = 14, height = 7.4, dpi = 300)
  ggsave(file.path(fig_dir, "Fig_S3F_strict_sensitive_within_genus_pairwise_species_heatmap.pdf"), p_pair, width = 14, height = 7.4)
}

summary_lines <- c(
  "# eCF presence/absence taxonomic visualization",
  "",
  "Input: strict-sensitive MAG eCF copy-number matrix with phenotypes.",
  "",
  "Presence definition: copy_number > 0. Copy number is used only as an auxiliary color scale for mean copy among eCF-positive MAGs.",
  "",
  "Outputs:",
  "- Fig_S3C_strict_sensitive_genus_eCF_presence_copy_deployment.png and .pdf",
  "- Fig_S3D_strict_sensitive_within_genus_species_eCF_deployment.png and .pdf",
  "- Fig_S3E_strict_sensitive_eCF_prevalence_range_rank.png and .pdf",
  "- Fig_S3F_strict_sensitive_within_genus_pairwise_species_heatmap.png and .pdf, when pairwise tests are available",
  "",
  "Statistical interpretation:",
  "- genus-level range is max genus prevalence minus min genus prevalence for each eCF;",
  "- species-level range is calculated separately within each genus for each genus-eCF test;",
  "- species pairwise tests are only drawn for eCFs with within-genus global FDR < 0.05."
)
writeLines(summary_lines, file.path(fig_dir, "CF_PRESENCE_VISUALIZATION_SUMMARY.md"))

message("eCF presence/absence figures written: ", fig_dir)

