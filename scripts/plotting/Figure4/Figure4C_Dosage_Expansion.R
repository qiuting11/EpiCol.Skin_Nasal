# ==============================================================================
# Figure4C_Dosage_Expansion.R
# ==============================================================================
# Positioning: Figure 4 - Dosage-focused expansion analysis at Broader Scale.
# Logic: Aligned with Figure 3B genera for consistent storytelling.
# Visualization: Size = Mean Copy (Pos), Color = Binned Mean Copy (Pos).
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(vroom)
  library(ggplot2)
  library(patchwork)
  library(scales)
})

script_dir <- "scripts/plotting/Figure3"
if (file.exists(file.path(script_dir, "Figure3_color_config.R"))) {
  source(file.path(script_dir, "Figure3_color_config.R"))
}

paths <- list(
  modules = "data/figure4/top0.2_I3.0_modules.txt",
  copy_counts = "data/external/true_cf_genome_copy_counts.csv",
  taxonomy = "data/external/bac120_taxonomy.tsv",
  f3b_deployment = "data/figure3/Figure3B_Genus_AllCF_Deployment.tsv",
  out_dir = "outputs/Figure4"
)
dir.create(paths$out_dir, recursive = TRUE, showWarnings = FALSE)

message("Syncing genera with Figure 3B...")
f3b_meta <- vroom(paths$f3b_deployment, show_col_types = FALSE) %>%
  distinct(Genus, HabitatGroup, n_genomes) %>%
  arrange(
    factor(HabitatGroup, levels = c("Skin-associated", "Shared", "Nasal-associated", "Other")),
    desc(n_genomes),
    Genus
  )

target_genera_clean <- f3b_meta$Genus
genus_levels <- rev(target_genera_clean)

message("Reading module info...")
modules_df <- read_delim(paths$modules, delim = "\t", show_col_types = FALSE) %>%
  mutate(
    module = factor(module, levels = module_order)
  ) %>%
  arrange(module, cf_id) %>%
  mutate(
    functional_theme = unname(functional_theme_map[as.character(module)]),
    functional_theme_short = unname(functional_theme_short[functional_theme])
  )

message("Reading taxonomy and mapping genera...")
tax_df <- vroom(paths$taxonomy, delim = "\t", col_names = c("genome_id", "taxonomy"), show_col_types = FALSE) %>%
  mutate(Genus = str_extract(taxonomy, "(?<=g__)[^;]+")) %>%
  filter(!is.na(Genus))

genus_totals <- tax_df %>%
  group_by(Genus) %>%
  summarise(total_genomes = n(), .groups = "drop")

target_genome_map <- tax_df %>%
  filter(Genus %in% target_genera_clean) %>%
  select(genome_id, Genus)

message("Reading DEDUPLICATED copy count data...")
copy_counts <- vroom(paths$copy_counts, show_col_types = FALSE)

message("Calculating mean copy number...")
plot_data_raw <- copy_counts %>%
  inner_join(target_genome_map, by = "genome_id")

metrics_df <- plot_data_raw %>%
  group_by(Genus, CF_id) %>%
  summarise(
    positive_genomes = n_distinct(genome_id),
    mean_copy_pos = mean(true_copy_number),
    .groups = "drop"
  )

missing_cfs_known <- c("CF15_19", "CF17_23", "CF9_12")
valid_cfs <- unique(metrics_df$CF_id)
cf_order <- modules_df %>% 
  filter(!cf_id %in% missing_cfs_known) %>%
  filter(cf_id %in% valid_cfs) %>%
  pull(cf_id) %>% 
  unique()

modules_df <- modules_df %>% filter(cf_id %in% cf_order)

theme_short_labels <- c(
  "Nutrient/substrate use" = "Nutrient",
  "Nucleotide/genome maintenance" = "Nucleotide",
  "Stress/envelope homeostasis" = "Stress/env.",
  "Translation/RNA regulation" = "Translation/RNA",
  "Communication signal" = "Signal",
  "Poorly characterized" = "Poorly char."
)

copy_bin_levels <- c("1", ">1-2", ">2-5", ">5-10", ">10-20", ">20-40", ">40")
copy_bin_colors <- c(
  "1"       = "#D9D9D9",
  ">1-2"    = "#FEE5D9",
  ">2-5"    = "#FCBBA1",
  ">5-10"   = "#FC9272",
  ">10-20"  = "#FB6A4A",
  ">20-40"  = "#DE2D26",
  ">40"     = "#A50F15"
)

plot_df <- metrics_df %>%
  filter(CF_id %in% cf_order) %>%
  tidyr::complete(Genus, CF_id = cf_order) %>%
  left_join(modules_df %>% select(CF_id = cf_id, module, functional_theme_short), by = "CF_id") %>%
  mutate(
    Genus = factor(Genus, levels = genus_levels),
    CF_id = factor(CF_id, levels = cf_order),
    functional_theme_plot = factor(unname(theme_short_labels[as.character(functional_theme_short)]), 
                                   levels = unname(theme_short_labels)),
    copy_bin = factor(case_when(
      mean_copy_pos > 0 & mean_copy_pos <= 1 ~ "1",
      mean_copy_pos > 1 & mean_copy_pos <= 2 ~ ">1-2",
      mean_copy_pos > 2 & mean_copy_pos <= 5 ~ ">2-5",
      mean_copy_pos > 5 & mean_copy_pos <= 10 ~ ">5-10",
      mean_copy_pos > 10 & mean_copy_pos <= 20 ~ ">10-20",
      mean_copy_pos > 20 & mean_copy_pos <= 40 ~ ">20-40",
      mean_copy_pos > 40 ~ ">40",
      TRUE ~ NA_character_
    ), levels = copy_bin_levels)
  )

# 6. Build panels
annotation_df <- modules_df %>%
  distinct(cf_id, module, functional_theme_short) %>%
  mutate(
    CF_id = factor(cf_id, levels = cf_order),
    functional_theme_plot = factor(unname(theme_short_labels[as.character(functional_theme_short)]), 
                                   levels = unname(theme_short_labels)),
    fill_color = unname(full_module_colors[as.character(module)]),
    label = str_replace(module, "module", "eCM")
  )

ann <- ggplot(annotation_df, aes(x = CF_id, y = "eCM")) +
  geom_tile(aes(fill = fill_color), width = 0.96, height = 0.82, color = "white", linewidth = 0.25) +
  geom_text(aes(label = label), size = 1.9, color = "white") +
  facet_grid(. ~ functional_theme_plot, scales = "free_x", space = "free_x") +
  scale_fill_identity() +
  labs(x = NULL, y = NULL, title = "Figure 4C | eCF Dosage (GTDB Scale) Aligned with Fig 3B Genera") +
  theme_minimal(base_size = 8) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(color = "black", size = 7),
    strip.text.x = element_text(size = 7, face = "bold"),
    strip.background = element_rect(fill = "grey96", color = "white"),
    plot.title = element_text(size = 12, face = "bold"),
    plot.margin = margin(5, 5, 0, 5)
  )

bubble <- ggplot(plot_df, aes(x = CF_id, y = Genus)) +
  geom_point(
    data = plot_df %>% filter(!is.na(copy_bin)),
    aes(color = copy_bin),
    size = 3.5,
    alpha = 0.94
  ) +
  facet_grid(. ~ functional_theme_plot, scales = "free_x", space = "free_x") +
  scale_x_discrete(labels = function(x) sub("^CF", "eCF", x)) +
  scale_color_manual(
    name = "Copy bin",
    values = copy_bin_colors,
    limits = copy_bin_levels,
    drop = FALSE
  ) +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 8) +
  theme(
    panel.grid.major = element_line(color = "grey92", linewidth = 0.25),
    panel.grid.minor = element_blank(),
    axis.text.y = element_text(color = "black", size = 7),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 7),
    strip.text.x = element_blank(),
    strip.background = element_blank(),
    plot.margin = margin(2, 5, 5, 5),
    legend.position = "bottom"
  )

summary_stats <- genus_totals %>%
  filter(Genus %in% target_genera_clean) %>%
  mutate(Genus = factor(Genus, levels = genus_levels))

cf_richness <- plot_df %>%
  group_by(Genus) %>%
  summarise(eCF_richness = n_distinct(CF_id[!is.na(copy_bin)]), .groups = "drop") %>%
  mutate(Genus = factor(Genus, levels = genus_levels))

track_theme <- theme_minimal(base_size = 8) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_text(size = 6),
    plot.margin = margin(2, 3, 5, 3)
  )

p_n_genomes <- ggplot(summary_stats, aes(x = total_genomes, y = Genus)) +
  geom_col(fill = "#8C8C8C", width = 0.72) +
  scale_x_log10(labels = label_number(scale_cut = cut_short_scale())) +
  scale_y_discrete(drop = FALSE) +
  labs(x = "GTDB size", y = NULL) +
  track_theme

p_richness <- ggplot(cf_richness, aes(x = eCF_richness, y = Genus)) +
  geom_col(fill = "#4C78A8", width = 0.72) +
  scale_y_discrete(drop = FALSE) +
  labs(x = "eCF detected", y = NULL) +
  track_theme

make_track_header <- function(title) {
  ggplot() +
    xlim(0, 1) + ylim(0, 1) +
    annotate("text", x = 0.5, y = 0.5, label = title, fontface = "bold", size = 3) +
    theme_void() +
    theme(plot.margin = margin(5, 3, 0, 3))
}

message("Assembling plot panels...")
body_height <- max(6.0, length(target_genera_clean) * 0.25)
left_panel <- ann / bubble + plot_layout(heights = c(0.55, body_height))
track_panel_1 <- make_track_header("GTDB Size") / p_n_genomes + plot_layout(heights = c(0.55, body_height))
track_panel_2 <- make_track_header("eCF Count") / p_richness + plot_layout(heights = c(0.55, body_height))

combined <- wrap_plots(left_panel, track_panel_1, track_panel_2, nrow = 1, widths = c(10, 0.8, 0.8), guides = "collect") &
  theme(legend.position = "bottom")

plot_height <- max(7.4, length(target_genera_clean) * 0.28)
ggsave(file.path(paths$out_dir, "Figure4C_Dosage_Expansion_Aligned.png"), combined, width = 16, height = plot_height, dpi = 300)
ggsave(file.path(paths$out_dir, "Figure4C_Dosage_Expansion_Aligned.pdf"), combined, width = 16, height = plot_height)

message("Figure4C Bubble Plot completed.")
