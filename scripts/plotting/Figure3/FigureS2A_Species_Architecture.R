suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
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


all_species_df <- read.delim(file.path(in_dir, "Figure3C_Species_AllCF_Architecture.tsv"),
                             sep = "\t", header = TRUE, check.names = FALSE)

literature_priority_species <- tibble::tribble(
  ~Species, ~selection_group, ~selection_reason, ~literature_support,
  "Staphylococcus epidermidis", "Skin/Nasal Staphylococcus", "core epithelial-associated skin and nasal Staphylococcus; retained from v3", "Liu2026_nasal; Liu2025_skin; CMR2026_review",
  "Staphylococcus aureus", "Skin/Nasal Staphylococcus", "clinically important skin/nasal carriage species; retained from v3", "Liu2026_nasal; CMR2026_review",
  "Staphylococcus hominis", "Skin Staphylococcus", "common epithelial-associated skin coagulase-negative Staphylococcus; retained from v3", "Liu2025_skin; CMR2026_review",
  "Staphylococcus capitis", "Skin Staphylococcus", "common epithelial-associated skin coagulase-negative Staphylococcus; retained from v3", "Liu2025_skin; CMR2026_review",
  "Staphylococcus warneri", "Skin Staphylococcus", "skin-associated coagulase-negative Staphylococcus; retained from v3", "Liu2025_skin; CMR2026_review",
  "Staphylococcus saprophyticus", "Skin Staphylococcus", "skin-side Staphylococcus retained for continuity from v3", "Liu2025_skin",
  "Staphylococcus simulans", "Skin Staphylococcus", "skin-side Staphylococcus retained for continuity from v3", "Liu2025_skin",
  "Staphylococcus auricularis", "Skin Staphylococcus", "skin-side Staphylococcus retained for continuity from v3", "Liu2025_skin",
  "Staphylococcus haemolyticus", "Skin Staphylococcus", "literature-supported epithelial-associated skin Staphylococcus; added in v4 despite n<5", "Liu2025_skin; CMR2026_review",
  "Staphylococcus lugdunensis", "Skin Staphylococcus", "clinically relevant epithelial-associated skin Staphylococcus; added in v4 despite n<5", "Liu2025_skin; CMR2026_review",
  "Staphylococcus saccharolyticus", "Skin Staphylococcus", "anaerobic skin-associated Staphylococcus; added in v4", "Liu2025_skin; CMR2026_review",
  "Staphylococcus pettenkoferi", "Skin Staphylococcus", "epithelial-associated skin Staphylococcus represented in current data; added in v4", "Liu2025_skin; CMR2026_review",
  "Staphylococcus cohnii", "Skin Staphylococcus", "skin-associated Staphylococcus represented in current data; added in v4", "Liu2025_skin; CMR2026_review",
  "Cutibacterium acnes", "Skin Cutibacterium", "canonical lipid-rich skin Cutibacterium; retained from v3", "Liu2025_skin; CMR2026_review",
  "Cutibacterium granulosum", "Skin Cutibacterium", "epithelial-associated skin Cutibacterium representative; retained from v3", "Liu2025_skin; CMR2026_review",
  "Cutibacterium modestum", "Skin Cutibacterium", "skin-side Cutibacterium represented in current data; retained from v3", "Liu2025_skin",
  "Cutibacterium namnetense", "Skin Cutibacterium", "skin-side Cutibacterium represented in current data; retained from v3", "Liu2025_skin",
  "Cutibacterium avidum", "Skin Cutibacterium", "epithelial-associated skin Cutibacterium; added in v4 despite n<5", "Liu2025_skin; CMR2026_review",
  "Corynebacterium accolens", "Nasal/shared Corynebacterium", "nasal-associated Corynebacterium; retained from v3", "Liu2026_nasal; CMR2026_review",
  "Corynebacterium propinquum", "Nasal/shared Corynebacterium", "nasal-enriched Corynebacterium; retained from v3", "Liu2026_nasal; CMR2026_review",
  "Corynebacterium pseudodiphtheriticum", "Nasal/shared Corynebacterium", "nasal-associated Corynebacterium; retained from v3", "Liu2026_nasal; CMR2026_review",
  "Corynebacterium striatum", "Skin/Nasal Corynebacterium", "epithelial-associated skin/nasal Corynebacterium; added in v4 despite n<5", "Liu2025_skin; Liu2026_nasal; CMR2026_review",
  "Corynebacterium kroppenstedtii", "Skin Corynebacterium", "lipophilic skin-associated Corynebacterium; added in v4 despite n<5", "Liu2025_skin; CMR2026_review",
  "Dolosigranulum pigrum", "Nasal commensal", "nasal commensal representative; retained from v3", "Liu2026_nasal",
  "Moraxella catarrhalis", "Nasal Moraxella", "nasal-associated Moraxella; retained from v3", "Liu2026_nasal",
  "Moraxella nonliquefaciens", "Nasal Moraxella", "nasal-associated Moraxella; retained from v3", "Liu2026_nasal",
  "Moraxella_A osloensis", "Skin/Nasal Moraxella", "Moraxella osloensis lineage represented in current data; added in v4", "Liu2026_nasal; CMR2026_review",
  "Rothia mucilaginosa", "Upper-airway/shared", "upper-airway-associated species represented in current data; added in v4", "Liu2026_nasal",
  "Streptococcus pneumoniae", "Nasal/airway Streptococcus", "nasal carriage pathogen/commensal context; added in v4 despite n<5", "Liu2026_nasal; CMR2026_review"
)

species_metrics <- all_species_df %>%
  distinct(Species, Genus, n_genomes, n_skin, n_nasal)

selected_species_metadata <- literature_priority_species %>%
  left_join(species_metrics, by = "Species") %>%
  mutate(
    available_in_current_data = !is.na(n_genomes),
    manual_low_n_keep = available_in_current_data & n_genomes < 5,
    selected_for_figure3c_v4 = available_in_current_data
  )

selected_species <- selected_species_metadata %>%
  filter(selected_for_figure3c_v4) %>%
  pull(Species)

df <- all_species_df %>%
  filter(Species %in% selected_species) %>%
  mutate(Species = factor(Species, levels = selected_species)) %>%
  arrange(Species, module, CF)

write.table(selected_species_metadata,
            file.path(out_dir, "Figure3C_KeySpecies_Rationale_v4.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(data.frame(SelectedSpecies = selected_species),
            file.path(out_dir, "Figure3C_KeySpecies_Selected_v4.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(df,
            file.path(out_dir, "Figure3C_KeySpecies_AllCF_Architecture_v4.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

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

# Key fix: Only keep CFs that exist in the current species list and define their order to remove empty columns.
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
  tidyr::complete(Species, CF = cf_order) %>%
  mutate(
    Species = factor(Species, levels = rev(selected_species)),
    CF = factor(as.character(CF), levels = cf_order),
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

species_summary <- df %>%
  distinct(Species, Genus, n_genomes, n_skin, n_nasal) %>%
  mutate(Species = factor(Species, levels = rev(selected_species)))

species_levels <- levels(species_summary$Species)

stack_df <- species_summary %>%
  select(Species, Skin = n_skin, Nasal = n_nasal) %>%
  pivot_longer(cols = c(Skin, Nasal), names_to = "Habitat", values_to = "n") %>%
  mutate(Habitat = factor(Habitat, levels = c("Nasal", "Skin")))

ann <- ggplot(annotation_df, aes(x = CF, y = annotation)) +
  geom_tile(aes(fill = fill_color), width = 0.96, height = 0.82, color = "white", linewidth = 0.25) +
  geom_text(aes(label = label), size = 1.9, color = "black") +
  facet_grid(. ~ functional_theme_plot, scales = "free_x", space = "free_x") +
  scale_fill_identity() +
  labs(x = NULL, y = NULL, title = "Figure 3C v4 | Literature-guided eCF architecture of representative species") +
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

bubble <- ggplot(plot_df, aes(x = CF, y = Species)) +
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

stack_track <- ggplot(stack_df, aes(x = n, y = Species, fill = Habitat)) +
  geom_col(width = 0.72) +
  scale_fill_manual(values = c(Skin = "#D9A441", Nasal = "#5B8DB8"), breaks = c("Skin", "Nasal"), name = "Genome source") +
  scale_y_discrete(drop = FALSE, limits = species_levels) +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 8) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_text(size = 6),
    plot.margin = margin(2, 3, 5, 3)
  )

track_header <- ggplot() +
  annotate("text", x = 0.5, y = 0.5, label = "Skin/Nasal\nn genomes", fontface = "bold", size = 3) +
  xlim(0, 1) +
  ylim(0, 1) +
  theme_void() +
  theme(plot.margin = margin(5, 3, 0, 3))

body_height <- max(6.4, length(selected_species) * 0.27)
left_panel <- ann / bubble + plot_layout(heights = c(0.55, body_height))
track_panel <- track_header / stack_track + plot_layout(heights = c(0.55, body_height))

combined <- wrap_plots(left_panel, track_panel, nrow = 1, widths = c(10, 1.15), guides = "collect") &
  theme(legend.position = "bottom")

plot_height <- max(7.4, length(selected_species) * 0.34)
ggsave(file.path(out_dir, "Figure3C_KeySpecies_AllCF_CMAnnotated_v4.png"),
       combined, width = 16, height = plot_height, dpi = 300)
ggsave(file.path(out_dir, "Figure3C_KeySpecies_AllCF_CMAnnotated_v4.pdf"),
       combined, width = 16, height = plot_height)

message("Figure3C v4 completed. Selected species: ", length(selected_species),
        ". Copy color uses fixed bins: ", paste(copy_bin_levels, collapse = ", "))
