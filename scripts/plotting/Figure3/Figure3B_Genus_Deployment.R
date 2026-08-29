suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(stringr)
})

script_dir <- "scripts/plotting/Figure3"
source(file.path(script_dir, "Figure3_color_config.R"))

in_dir <- "data/figure3"
out_dir <- "outputs/Figure3"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)


df <- read.delim(file.path(in_dir, "Figure3B_Genus_AllCF_Deployment.tsv"),
                 sep = "\t", header = TRUE, check.names = FALSE)

theme_order <- names(functional_theme_colors)
theme_short_levels <- unname(functional_theme_short[theme_order])
theme_short_labels <- c(
  "Nutrient/substrate use" = "Nutrient",
  "Nucleotide/genome maintenance" = "Nucleotide",
  "Stress/envelope homeostasis" = "Stress/env.",
  "Translation/RNA regulation" = "Translation/RNA",
  "Communication signal" = "Signal",
  "Poorly characterized" = "Poorly char."
)

genus_summary <- df %>%
  distinct(Genus, HabitatGroup, n_genomes, n_species, n_skin, n_nasal, CF_richness,
           mean_total_CF_copy_per_genome, main_selection_reason) %>%
  arrange(
    factor(HabitatGroup, levels = c("Skin-associated", "Shared", "Nasal-associated", "Other")),
    desc(n_genomes),
    Genus
  ) %>%
  mutate(Genus = as.character(Genus))

genus_levels <- rev(genus_summary$Genus)

copy_bin_levels <- c("1", ">1-1.5", ">1.5-2", ">2-3", ">3-5", ">5-10", ">10")
copy_bin_colors <- c(
  "1"       = "#D9D9D9",
  ">1-1.5"  = "#FEE0D2",
  ">1.5-2"  = "#FCBBA1",
  ">2-3"    = "#FC9272",
  ">3-5"    = "#FB6A4A",
  ">5-10"   = "#DE2D26",
  ">10"     = "#A50F15"
)

# Key fix: Only keep CFs that exist in the current genus list and define their order to remove empty columns.
valid_cfs <- unique(df$CF[df$mean_copy_pos > 0])
cf_meta <- df %>%
  filter(CF %in% valid_cfs) %>%
  distinct(CF, module, functional_theme, functional_theme_short) %>%
  arrange(factor(as.character(module), levels = module_order), factor(as.character(CF), levels = unique(as.character(CF)))) %>%
  mutate(
    CF = factor(as.character(CF), levels = unique(as.character(CF))),
    module = as.character(module),
    functional_theme_short = factor(functional_theme_short, levels = theme_short_levels),
    functional_theme_plot = factor(unname(theme_short_labels[as.character(functional_theme_short)]),
                                   levels = unname(theme_short_labels[theme_short_levels])),
    module_fill = unname(full_module_colors[module])
  )

cf_order <- levels(cf_meta$CF)

annotation_df <- cf_meta %>%
  transmute(
    CF,
    functional_theme_plot,
    annotation = factor("eCM", levels = "eCM"),
    fill_color = module_fill,
    label = str_replace(module, "module", "eCM")
  )

plot_df <- df %>%
  filter(CF %in% cf_order) %>%
  tidyr::complete(Genus, CF = cf_order) %>%
  left_join(genus_summary %>% select(Genus, HabitatGroup), by = "Genus") %>% # Information lost by 'complete' needs to be recovered.
  mutate(
    Genus = factor(as.character(Genus), levels = genus_levels),
    CF = factor(as.character(CF), levels = cf_order),
    module = factor(as.character(module), levels = module_order),
    functional_theme_short = factor(functional_theme_short, levels = theme_short_levels),
    functional_theme_plot = factor(unname(theme_short_labels[as.character(functional_theme_short)]),
                                   levels = unname(theme_short_labels[theme_short_levels])),
    copy_bin = factor(case_when(
      mean_copy_pos > 0 & mean_copy_pos <= 1 ~ "1",
      mean_copy_pos > 1 & mean_copy_pos <= 1.5 ~ ">1-1.5",
      mean_copy_pos > 1.5 & mean_copy_pos <= 2 ~ ">1.5-2",
      mean_copy_pos > 2 & mean_copy_pos <= 3 ~ ">2-3",
      mean_copy_pos > 3 & mean_copy_pos <= 5 ~ ">3-5",
      mean_copy_pos > 5 & mean_copy_pos <= 10 ~ ">5-10",
      mean_copy_pos > 10 ~ ">10",
      TRUE ~ NA_character_
    ), levels = copy_bin_levels)
  )

ann <- ggplot(annotation_df, aes(x = CF, y = annotation)) +
  geom_tile(aes(fill = fill_color), width = 0.96, height = 0.82, color = "white", linewidth = 0.25) +
  geom_text(aes(label = label), size = 1.9, color = "black") +
  facet_grid(. ~ functional_theme_plot, scales = "free_x", space = "free_x") +
  scale_fill_identity() +
  labs(x = NULL, y = NULL, title = "Figure 3B | Genus-level all-eCF deployment") +
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

bubble <- ggplot(plot_df, aes(x = CF, y = Genus)) +
  geom_point(
    data = plot_df %>% filter(mean_copy_pos > 0, prevalence > 0),
    aes(size = prevalence, color = copy_bin),
    alpha = 0.94
  ) +
  facet_grid(. ~ functional_theme_plot, scales = "free_x", space = "free_x") +
  scale_size_continuous(name = "Prevalence", range = c(0.2, 4.2), limits = c(0, 1), breaks = c(0.25, 0.5, 0.75, 1)) +
  scale_x_discrete(labels = function(x) sub("^CF", "eCF", x)) +
  scale_color_manual(
    name = "Mean copy\n(positive)",
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

summary_df <- genus_summary %>%
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

make_track <- function(data, xvar, fill) {
  ggplot(data, aes(x = .data[[xvar]], y = Genus)) +
    geom_col(fill = fill, width = 0.72) +
    scale_y_discrete(drop = FALSE, limits = genus_levels) +
    labs(x = NULL, y = NULL) +
    track_theme
}

make_track_header <- function(title) {
  ggplot() +
    annotate("text", x = 0.5, y = 0.5, label = title, fontface = "bold", size = 3) +
    xlim(0, 1) +
    ylim(0, 1) +
    theme_void() +
    theme(plot.margin = margin(5, 3, 0, 3))
}

p_n_genomes <- make_track(summary_df, "n_genomes", "#8C8C8C")
p_richness <- make_track(summary_df, "CF_richness", "#4C78A8")

body_height <- max(6.0, length(genus_levels) * 0.25)
left_panel <- ann / bubble + plot_layout(heights = c(0.55, body_height))
track_panel_1 <- make_track_header("Genomes") / p_n_genomes + plot_layout(heights = c(0.55, body_height))
track_panel_2 <- make_track_header("eCF rich.") / p_richness + plot_layout(heights = c(0.55, body_height))

combined <- wrap_plots(left_panel, track_panel_1, track_panel_2, nrow = 1, widths = c(10, 0.8, 0.8), guides = "collect") &
  theme(legend.position = "bottom")

plot_height <- max(7.4, length(genus_levels) * 0.28)
ggsave(file.path(out_dir, "Figure3B_Genus_AllCF_Bubble_CMOnly_v3.png"),
       combined, width = 16, height = plot_height, dpi = 300)
ggsave(file.path(out_dir, "Figure3B_Genus_AllCF_Bubble_CMOnly_v3.pdf"),
       combined, width = 16, height = plot_height)

message("Figure3B v3 completed with eCF/eCM display labels. Copy color uses fixed bins: ", paste(copy_bin_levels, collapse = ", "))
