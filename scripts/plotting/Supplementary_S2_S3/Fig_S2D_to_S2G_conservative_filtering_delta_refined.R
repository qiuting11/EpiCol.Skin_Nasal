#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(readr)
  library(stringr)
  library(tidyr)
})

repo_root <- Sys.getenv("PAPER_SUBMISSION_CODE", unset = normalizePath(getwd(), winslash = "/", mustWork = TRUE))
root <- dirname(repo_root)
old_dir <- file.path(root, "R", "gephe_R", "results", "figures", "Figure3_Final")
new_dir <- file.path(root, "R", "gephe_R", "results", "figures", "Figure3_human_pre11_paper_id50_k1_original_style")
module_path <- file.path(repo_root, "data", "figure3", "top0.2_I3.0_modules.txt")
color_config <- file.path(repo_root, "scripts", "plotting", "Figure3", "Figure3_color_config.R")
out_dir <- file.path(repo_root, "outputs", "Supplementary_S2_S3")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
plot_table_dir <- file.path(out_dir, "plot_tables")
dir.create(plot_table_dir, recursive = TRUE, showWarnings = FALSE)
source(color_config)

mods <- read_tsv(module_path, show_col_types = FALSE) %>%
  transmute(
    CF = cf_id,
    eCM = module,
    eCM_num = as.integer(str_remove(module, "module")),
    cf_major = as.integer(str_match(cf_id, "^CF([0-9]+)")[, 2]),
    cf_minor = suppressWarnings(as.integer(str_match(cf_id, "_([0-9]+)$")[, 2])),
    cf_minor = if_else(is.na(cf_minor), -1L, cf_minor),
    eCF_label = str_replace(cf_id, "^CF", "eCF"),
    eCM_label = str_replace(module, "^module", "eCM"),
    module_fill = unname(full_module_colors[module])
  ) %>%
  arrange(eCM_num, cf_major, cf_minor, CF) %>%
  mutate(x_pos = row_number())

cf_order <- mods$CF

read_deployment <- function(path, row_col) {
  read_tsv(path, show_col_types = FALSE) %>%
    select(all_of(row_col), CF, prevalence, mean_copy_all)
}

build_delta <- function(old_path, new_path, row_col, selection_path, selection_col, out_prefix) {
  old <- read_deployment(old_path, row_col)
  new <- read_deployment(new_path, row_col)
  row_order <- read_tsv(selection_path, show_col_types = FALSE) %>%
    pull(all_of(selection_col)) %>%
    as.character() %>%
    unique()

  delta <- old %>%
    inner_join(new, by = c(row_col, "CF"), suffix = c("_broad", "_strict")) %>%
    mutate(
      delta_prevalence = prevalence_strict - prevalence_broad,
      delta_mean_copy = mean_copy_all_strict - mean_copy_all_broad,
      row_label = .data[[row_col]]
    ) %>%
    filter(CF %in% cf_order, row_label %in% row_order) %>%
    mutate(
      CF = factor(CF, levels = cf_order),
      row_label = factor(row_label, levels = rev(row_order))
    )

  write_tsv(delta, file.path(plot_table_dir, paste0(out_prefix, ".tsv")))
  delta
}

make_ecm_strip <- function() {
  ggplot(mods, aes(x = x_pos, y = 1, fill = module_fill)) +
    geom_tile(width = 0.96, height = 0.62, color = "white", linewidth = 0.20) +
    geom_text(aes(label = eCM_label), size = 1.6, color = "black") +
    scale_fill_identity() +
    scale_x_continuous(breaks = mods$x_pos, labels = NULL, expand = expansion(mult = c(0.004, 0.004))) +
    scale_y_continuous(expand = c(0, 0)) +
    theme_void(base_size = 8.2, base_family = "sans") +
    theme(plot.margin = margin(0, 16, 0, 44))
}

heatmap_plot <- function(delta, title, file_base, height) {
  plot_df <- delta %>% left_join(mods %>% select(CF, x_pos, eCF_label), by = "CF")
  p_title <- ggplot() +
    annotate("text", x = 0, y = 0.54, label = title, hjust = 0, fontface = "bold", size = 3.0) +    xlim(0, 1) + ylim(0, 1) + theme_void() +
    theme(plot.margin = margin(2, 18, 0, 56))
  p <- ggplot(plot_df, aes(x = x_pos, y = row_label, fill = delta_prevalence)) +
    geom_tile(color = "white", linewidth = 0.12) +
    scale_x_continuous(breaks = mods$x_pos, labels = mods$eCF_label, expand = expansion(mult = c(0.004, 0.004))) +
    scale_fill_gradient2(
      low = "#3B6FB6", mid = "#F8F8F8", high = "#B64A4A", midpoint = 0,
      limits = c(-1, 1), name = "Delta prevalence", labels = scales::percent_format(accuracy = 1)
    ) +
    labs(x = NULL, y = NULL) +
    theme_classic(base_size = 8.2, base_family = "sans") +
    theme(
      panel.grid = element_blank(),
      axis.line = element_line(linewidth = 0.25, color = "black"),
      axis.ticks = element_line(linewidth = 0.25, color = "black"),
      axis.ticks.length = unit(1.2, "mm"),
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6.2, color = "black"),
      axis.text.y = element_text(size = 6.8, color = "black"),
      legend.position = "right",
      legend.title = element_text(size = 7.2),
      legend.text = element_text(size = 6.6),
      plot.margin = margin(0, 18, 5, 5)
    )
  combined <- p_title / make_ecm_strip() / p + plot_layout(heights = c(0.10, 0.065, 1))
  ggsave(file.path(out_dir, paste0(file_base, ".png")), combined, width = 15.4, height = height, dpi = 300)
  ggsave(file.path(out_dir, paste0(file_base, ".pdf")), combined, width = 15.4, height = height)
}

summary_plot <- function(delta, title, file_base) {
  summary_df <- delta %>%
    group_by(CF) %>%
    summarise(
      mean_delta_prevalence = mean(delta_prevalence, na.rm = TRUE),
      mean_delta_copy = mean(delta_mean_copy, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    left_join(mods %>% select(CF, x_pos, eCF_label), by = "CF") %>%
    arrange(x_pos)
  write_tsv(summary_df, file.path(plot_table_dir, paste0(file_base, ".tsv")))

  make_bar <- function(metric, y_label) {
    ggplot(summary_df, aes(x = x_pos, y = .data[[metric]], fill = .data[[metric]] >= 0)) +
      geom_col(width = 0.78) +
      geom_hline(yintercept = 0, linewidth = 0.25, color = "black") +
      scale_fill_manual(values = c(`TRUE` = "#B64A4A", `FALSE` = "#3B6FB6"), guide = "none") +
      scale_x_continuous(breaks = mods$x_pos, labels = NULL, expand = expansion(mult = c(0.004, 0.004))) +
      labs(x = NULL, y = y_label) +
      theme_classic(base_size = 8.2, base_family = "sans") +
      theme(
        panel.grid = element_blank(),
        axis.line = element_line(linewidth = 0.25, color = "black"),
        axis.ticks = element_line(linewidth = 0.25, color = "black"),
        axis.ticks.length = unit(1.2, "mm"),
        axis.text.y = element_text(size = 6.8, color = "black"),
        axis.title.y = element_text(size = 7.2, color = "black"),
        plot.margin = margin(1, 12, 1, 44)
      )
  }
  p_title <- ggplot() + annotate("text", x = 0, y = 0.5, label = title, hjust = 0, fontface = "bold", size = 3.0) +
    xlim(0, 1) + ylim(0, 1) + theme_void(base_family = "sans") + theme(plot.margin = margin(0, 12, 0, 44))
  p_x <- ggplot(summary_df, aes(x = x_pos, y = 1)) +
    scale_x_continuous(breaks = mods$x_pos, labels = mods$eCF_label, expand = expansion(mult = c(0.004, 0.004))) +
    scale_y_continuous(expand = c(0, 0)) +
    labs(x = NULL, y = NULL) + theme_void(base_size = 8.2, base_family = "sans") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6.2, color = "black"), plot.margin = margin(0, 12, 5, 56))
  combined <- p_title / make_ecm_strip() / make_bar("mean_delta_prevalence", "Mean delta
prevalence") / make_bar("mean_delta_copy", "Mean delta
copy") / p_x +
    plot_layout(heights = c(0.09, 0.065, 1, 1, 0.24))
  ggsave(file.path(out_dir, paste0(file_base, ".png")), combined, width = 13.6, height = 6.4, dpi = 300)
  ggsave(file.path(out_dir, paste0(file_base, ".pdf")), combined, width = 13.6, height = 6.4)
}

genus_delta <- build_delta(
  file.path(old_dir, "Figure3B_Genus_AllCF_Deployment.tsv"),
  file.path(new_dir, "Figure3B_Genus_AllCF_Deployment.tsv"),
  "Genus",
  file.path(new_dir, "Figure3B_SelectedGenera_Metadata.tsv"),
  "Genus",
  "Fig_S2D_S2E_delta_data"
)
heatmap_plot(genus_delta, "Genus-level prevalence change after conservative filtering", "Fig_S2D_genus_delta_prevalence_after_conservative_filtering", 8.2)
summary_plot(genus_delta, "Genus-level eCF summary after conservative filtering", "Fig_S2E_genus_delta_cf_summary_after_conservative_filtering")

species_delta <- build_delta(
  file.path(old_dir, "Figure3C_KeySpecies_AllCF_Architecture_v4.tsv"),
  file.path(new_dir, "Figure3C_KeySpecies_AllCF_Architecture_v4.tsv"),
  "Species",
  file.path(new_dir, "Figure3C_KeySpecies_Selected_v4.tsv"),
  "SelectedSpecies",
  "Fig_S2F_S2G_delta_data"
)
heatmap_plot(species_delta, "Species-level prevalence change after conservative filtering", "Fig_S2F_species_delta_prevalence_after_conservative_filtering", 9.4)
summary_plot(species_delta, "Species-level eCF summary after conservative filtering", "Fig_S2G_species_delta_cf_summary_after_conservative_filtering")

message("Fig. S2D-S2G files written to: ", out_dir)



