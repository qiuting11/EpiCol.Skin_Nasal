#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(cowplot)
  library(ggplot2)
})

repo_root <- Sys.getenv("PAPER_SUBMISSION_CODE", unset = normalizePath(getwd(), winslash = "/", mustWork = TRUE))
root <- file.path(repo_root, "data", "external", "supplementary_s2_s3", "human_hq_mag_cf_matrix")
tables_dir <- file.path(root, "controlled_parameter_comparison", "tables")
out_dir <- file.path(repo_root, "outputs", "Supplementary_S2_S3")
presence_tables_dir <- file.path(root, "main_strict_sensitive_result", "cf_presence_absence_analysis", "tables")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

source(file.path(repo_root, "scripts", "plotting", "Figure3", "Figure3_color_config.R"))

read_tsv_local <- function(path) {
  read.delim(path, sep = "\t", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
}

save_plot <- function(plot, path_base, width, height) {
  ggsave(paste0(path_base, ".pdf"), plot, width = width, height = height, units = "in", device = cairo_pdf)
  ggsave(paste0(path_base, ".png"), plot, width = width, height = height, units = "in", dpi = 300)
}

overall <- read_tsv_local(file.path(tables_dir, "controlled_three_way_cf_summary_long.tsv"))
cf_wide <- read_tsv_local(file.path(tables_dir, "controlled_three_way_cf_summary_wide.tsv"))
genus <- read_tsv_local(file.path(tables_dir, "controlled_three_way_genus_summary_long.tsv"))
genus_wide <- read_tsv_local(file.path(tables_dir, "controlled_three_way_genus_summary_wide.tsv"))

method_label_map <- c("strict-fast" = "Strict-fast", "strict-sensitive" = "Strict-sensitive", "paper-like" = "Broad homolog-recovery")
method_order <- c("Strict-fast", "Strict-sensitive", "Broad homolog-recovery")
cf_reference <- read_tsv_local(file.path(presence_tables_dir, "cf_presence_taxonomy_summary_for_plot.tsv"))
cf_reference <- cf_reference[!duplicated(cf_reference$CF_id), c("CF_id", "eCM", "eCM_num")]
cf_reference$cf_major <- as.integer(sub("^CF([0-9]+).*$", "\\1", cf_reference$CF_id))
cf_reference$cf_minor <- -1L
has_minor <- grepl("_", cf_reference$CF_id)
cf_reference$cf_minor[has_minor] <- as.integer(sub("^CF[0-9]+_([0-9]+)$", "\\1", cf_reference$CF_id[has_minor]))
cf_reference <- cf_reference[order(cf_reference$eCM_num, cf_reference$cf_major, cf_reference$cf_minor, cf_reference$CF_id), ]
cf_order <- cf_reference$CF_id
extra_cf <- setdiff(unique(overall$CF_id), cf_order)
if (length(extra_cf) > 0) {
  cf_order <- c(cf_order, sort(extra_cf))
}
cf_labels <- setNames(sub("^CF", "eCF", cf_order), cf_order)

module_bar <- cf_reference
module_bar$CF_id <- factor(module_bar$CF_id, levels = cf_order)
module_bar$eCM_label <- sub("^module", "eCM", module_bar$eCM)
module_bar$module_fill <- unname(full_module_colors[module_bar$eCM])
module_bar$module_fill[is.na(module_bar$module_fill)] <- "#D9D9D9"

overall$method <- method_label_map[overall$method]
overall$method[is.na(overall$method)] <- overall$method_id[is.na(overall$method)]
overall$method <- factor(overall$method, levels = method_order)
overall$CF_id <- factor(overall$CF_id, levels = cf_order)
genus$method <- method_label_map[genus$method]
genus$method[is.na(genus$method)] <- genus$method_id[is.na(genus$method)]
genus$method <- factor(genus$method, levels = method_order)
genus$CF_id <- factor(genus$CF_id, levels = cf_order)
genus_wide$CF_id <- factor(genus_wide$CF_id, levels = cf_order)

genus_label_map <- c(
  "Staphylococcus" = "Staphylococcus",
  "Corynebacterium" = "Corynebacterium",
  "Cutibacterium_or_Propionibacterium" = "Cutibacterium",
  "Cutibacterium" = "Cutibacterium"
)
genus_order <- c("Cutibacterium", "Corynebacterium", "Staphylococcus")
genus$target_genus_display <- genus_label_map[genus$target_genus_group]
genus$target_genus_display[is.na(genus$target_genus_display)] <- genus$target_genus_group[is.na(genus$target_genus_display)]
genus$target_genus_display <- factor(genus$target_genus_display, levels = genus_order)
genus_wide$target_genus_display <- genus_label_map[genus_wide$target_genus_group]
genus_wide$target_genus_display[is.na(genus_wide$target_genus_display)] <- genus_wide$target_genus_group[is.na(genus_wide$target_genus_display)]
genus_wide$target_genus_display <- factor(genus_wide$target_genus_display, levels = genus_order)

theme_cmp <- theme_bw(base_size = 9) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    strip.background = element_rect(fill = "grey92", color = "grey70"),
    legend.position = "right"
  )

theme_cmp_ecm <- theme_cmp +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7, color = "black"),
    axis.text.y = element_text(size = 8, color = "black"),
    strip.text.y = element_text(face = "bold", size = 8),
    plot.title = element_text(face = "bold", size = 11)
  )

make_module_strip <- function() {
  ggplot(module_bar, aes(x = CF_id, y = 1, fill = module_fill)) +
    geom_tile(color = "white", linewidth = 0.25) +
    geom_text(aes(label = eCM_label), size = 2, color = "black") +
    scale_fill_identity() +
    scale_x_discrete(limits = cf_order, labels = cf_labels, drop = FALSE, expand = expansion(add = 0.6)) +
    scale_y_continuous(expand = c(0, 0)) +
    labs(x = NULL, y = NULL) +
    theme_void(base_size = 8) +
    theme(plot.margin = margin(0, 50, 0, 72))
}

save_heatmap_with_module <- function(heat_plot, path_base, width, height, title = NULL, subtitle = NULL) {
  heat_plot <- heat_plot + labs(title = NULL, subtitle = NULL)
  if (!is.null(title)) {
    title_plot <- cowplot::ggdraw() +
      cowplot::draw_label(title, x = 0, y = ifelse(is.null(subtitle), 0.55, 0.68), hjust = 0, fontface = "bold", size = 11) +
      theme(plot.margin = margin(2, 48, 0, 72))
    if (!is.null(subtitle)) {
      title_plot <- title_plot + cowplot::draw_label(subtitle, x = 0, y = 0.25, hjust = 0, size = 8, color = "grey30")
    }
    p <- cowplot::plot_grid(title_plot, make_module_strip(), heat_plot, ncol = 1, rel_heights = c(0.16, 0.10, 1), align = "v", axis = "lr")
  } else {
    p <- cowplot::plot_grid(make_module_strip(), heat_plot, ncol = 1, rel_heights = c(0.10, 1), align = "v", axis = "lr")
  }
  save_plot(p, path_base, width, height)
}

p_overall_prev <- ggplot(overall, aes(x = CF_id, y = method, fill = prevalence)) +
  geom_tile(color = "white", linewidth = 0.25) +
  scale_fill_gradient(low = "white", high = "#1B7837", limits = c(0, 1), name = "Prevalence") +
  scale_x_discrete(limits = cf_order, labels = cf_labels, drop = FALSE) +
  labs(x = "eCF", y = NULL, title = "Overall eCF prevalence under controlled DIAMOND parameter sets") +
  theme_cmp

delta_prev <- rbind(
  data.frame(CF_id = cf_wide$CF_id, contrast = "Strict-sensitive minus Strict-fast", delta = cf_wide$delta_prevalence_sensitive_minus_fast),
  data.frame(CF_id = cf_wide$CF_id, contrast = "Broad homolog-recovery minus Strict-sensitive", delta = cf_wide$delta_prevalence_paper_minus_sensitive)
)
delta_prev$CF_id <- factor(delta_prev$CF_id, levels = cf_order)
delta_prev$contrast <- factor(delta_prev$contrast, levels = c("Strict-sensitive minus Strict-fast", "Broad homolog-recovery minus Strict-sensitive"))

p_delta_prev <- ggplot(delta_prev, aes(x = CF_id, y = contrast, fill = delta)) +
  geom_tile(color = "white", linewidth = 0.25) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0, limits = c(-1, 1), name = "Delta prevalence") +
  scale_x_discrete(limits = cf_order, labels = cf_labels, drop = FALSE) +
  labs(x = "eCF", y = NULL, title = "Controlled eCF prevalence changes") +
  theme_cmp

delta_copy <- rbind(
  data.frame(CF_id = cf_wide$CF_id, contrast = "Strict-sensitive minus Strict-fast", delta_total_copy = cf_wide$delta_total_copy_sensitive_minus_fast),
  data.frame(CF_id = cf_wide$CF_id, contrast = "Broad homolog-recovery minus Strict-sensitive", delta_total_copy = cf_wide$delta_total_copy_paper_minus_sensitive)
)
delta_copy$CF_id <- factor(delta_copy$CF_id, levels = cf_order)
delta_copy$contrast <- factor(delta_copy$contrast, levels = c("Strict-sensitive minus Strict-fast", "Broad homolog-recovery minus Strict-sensitive"))

p_delta_copy <- ggplot(delta_copy, aes(x = CF_id, y = delta_total_copy, fill = contrast)) +
  geom_col(position = "dodge", width = 0.75) +
  geom_hline(yintercept = 0, linewidth = 0.25) +
  labs(x = "eCF", y = "Delta total copy", title = "Controlled eCF total copy-number changes") +
  theme_cmp +
  theme(legend.position = "top")

p_genus_prev <- ggplot(genus, aes(x = CF_id, y = target_genus_display, fill = prevalence)) +
  geom_tile(color = "white", linewidth = 0.25) +
  facet_grid(method ~ ., switch = "y") +
  scale_fill_gradient(low = "white", high = "#1B7837", limits = c(0, 1), name = "Prevalence") +
  scale_x_discrete(limits = cf_order, labels = cf_labels, drop = FALSE, expand = expansion(add = 0.6)) +
  labs(x = "eCF", y = "Genus", title = "Genus-stratified eCF prevalence under controlled DIAMOND parameter sets") +
  theme_cmp_ecm

genus$log10_total_copy <- log10(genus$total_copy + 1)
p_genus_copy <- ggplot(genus, aes(x = CF_id, y = target_genus_display, fill = log10_total_copy)) +
  geom_tile(color = "white", linewidth = 0.25) +
  facet_grid(method ~ ., switch = "y") +
  scale_fill_gradient(low = "white", high = "#762A83", name = "log10(total copy + 1)") +
  scale_x_discrete(limits = cf_order, labels = cf_labels, drop = FALSE, expand = expansion(add = 0.6)) +
  labs(x = "eCF", y = "Genus", title = "Genus-stratified eCF total copy number under controlled DIAMOND parameter sets") +
  theme_cmp_ecm

genus_delta <- rbind(
  data.frame(target_genus_group = genus_wide$target_genus_group, CF_id = genus_wide$CF_id, contrast = "Strict-sensitive minus Strict-fast", delta = genus_wide$delta_prevalence_sensitive_minus_fast),
  data.frame(target_genus_group = genus_wide$target_genus_group, CF_id = genus_wide$CF_id, contrast = "Broad homolog-recovery minus Strict-sensitive", delta = genus_wide$delta_prevalence_paper_minus_sensitive)
)
genus_delta$CF_id <- factor(genus_delta$CF_id, levels = cf_order)
genus_delta$contrast <- factor(genus_delta$contrast, levels = c("Strict-sensitive minus Strict-fast", "Broad homolog-recovery minus Strict-sensitive"))
genus_delta$target_genus_display <- genus_label_map[genus_delta$target_genus_group]
genus_delta$target_genus_display[is.na(genus_delta$target_genus_display)] <- genus_delta$target_genus_group[is.na(genus_delta$target_genus_display)]
genus_delta$target_genus_display <- factor(genus_delta$target_genus_display, levels = genus_order)

p_genus_delta <- ggplot(genus_delta, aes(x = CF_id, y = target_genus_display, fill = delta)) +
  geom_tile(color = "white", linewidth = 0.25) +
  facet_grid(contrast ~ ., switch = "y") +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0, limits = c(-1, 1), name = "Delta prevalence") +
  scale_x_discrete(limits = cf_order, labels = cf_labels, drop = FALSE, expand = expansion(add = 0.6)) +
  labs(x = "eCF", y = "Genus", title = "Genus-stratified prevalence changes") +
  theme_cmp_ecm


copy_bin_colors <- c(
  "1" = "#D9D9D9",
  ">1-1.25" = "#FEE0D2",
  ">1.25-1.5" = "#FCBBA1",
  ">1.5-2" = "#FC9272",
  ">2-3" = "#FB6A4A",
  ">3" = "#A50F15"
)

genus$present_for_plot <- genus$prevalence > 0
genus$copy_bin <- NA_character_
genus$copy_bin[genus$present_for_plot & genus$mean_copy_present <= 1] <- "1"
genus$copy_bin[genus$present_for_plot & genus$mean_copy_present > 1 & genus$mean_copy_present <= 1.25] <- ">1-1.25"
genus$copy_bin[genus$present_for_plot & genus$mean_copy_present > 1.25 & genus$mean_copy_present <= 1.5] <- ">1.25-1.5"
genus$copy_bin[genus$present_for_plot & genus$mean_copy_present > 1.5 & genus$mean_copy_present <= 2] <- ">1.5-2"
genus$copy_bin[genus$present_for_plot & genus$mean_copy_present > 2 & genus$mean_copy_present <= 3] <- ">2-3"
genus$copy_bin[genus$present_for_plot & genus$mean_copy_present > 3] <- ">3"
genus$copy_bin <- factor(genus$copy_bin, levels = names(copy_bin_colors))

p_genus_bubble <- ggplot(genus[genus$present_for_plot, ], aes(x = CF_id, y = target_genus_display)) +
  geom_point(aes(size = prevalence, color = copy_bin), alpha = 0.95) +
  facet_grid(method ~ ., switch = "y") +
  scale_color_manual(
    name = "Mean copy\npresent MAGs",
    values = copy_bin_colors,
    limits = names(copy_bin_colors),
    drop = FALSE
  ) +
  scale_size_continuous(
    name = "eCF-positive MAG prevalence",
    range = c(0.6, 6.2),
    limits = c(0, 1),
    breaks = c(0.25, 0.5, 0.75, 1),
    labels = scales::percent_format(accuracy = 1)
  ) +
  scale_x_discrete(limits = cf_order, labels = cf_labels, drop = FALSE, expand = expansion(add = 0.6)) +
  labs(
    x = "eCF",
    y = "Genus",
    title = "Genus-stratified eCF prevalence and copy number under controlled DIAMOND parameter sets",
    subtitle = "Point size shows eCF-positive MAG prevalence; color shows mean copy number among eCF-positive MAGs."
  ) +
  theme_cmp_ecm +
  theme(
    panel.grid.major = element_line(color = "grey91", linewidth = 0.25),
    legend.position = "bottom",
    legend.box = "horizontal",
    plot.subtitle = element_text(size = 8, color = "grey30")
  )
save_heatmap_with_module(p_genus_bubble, file.path(out_dir, "Fig_S3A_controlled_genus_eCF_prevalence_copy_bubble"), 12.2, 7.0, title = "Genus-stratified eCF prevalence and copy number under controlled DIAMOND parameter sets", subtitle = "Point size shows eCF-positive MAG prevalence; color shows mean copy number among eCF-positive MAGs.")
writeLines(paste0("Fig. S3A file written: ", out_dir))
