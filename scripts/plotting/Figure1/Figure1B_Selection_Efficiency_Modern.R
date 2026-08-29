library(tidyverse)
library(grid)

# =========================
# 1. Source files
# =========================
all_human_file <- paste0(
  "data/figure1/",
  "all_human_mag_phenotypes.final_columns.tsv"
)
hq_report <- "data/figure1/00_merged_HQ_quality_report_with_rna.tsv"
final_metadata_file <- "data/figure1/metadata_for_analysis.tsv"

all_human <- read_tsv(all_human_file, show_col_types = FALSE)
hq_data <- read_tsv(hq_report, show_col_types = FALSE)
final_metadata <- read_tsv(final_metadata_file, show_col_types = FALSE)

# =========================
# 2. Epithelial-associated source classification
# =========================
classify_epithelial_source <- function(body_region, microenvironment) {
  case_when(
    body_region == "Nasal cavity" ~ "Nasal cavity",
    microenvironment == "dry" ~ "Dry skin",
    microenvironment == "sebaceous" ~ "Sebaceous skin",
    microenvironment == "moist" ~ "Moist skin",
    microenvironment == "toenail" ~ "Toenail",
    TRUE ~ "Unknown"
  )
}

make_epithelial_counts <- function(data, stage) {
  data %>%
    mutate(
      Stage = stage,
      Classification = "Epithelial-associated",
      Source = classify_epithelial_source(body_region, Microenvironments)
    ) %>%
    count(Stage, Classification, Source, name = "Value")
}

hq_epithelial_mags <- hq_data %>%
  filter(
    str_starts(check_Source, "UHSG") |
      str_starts(check_Source, "HSMG") |
      str_starts(check_Source, "nose")
  ) %>%
  distinct(MAG)

epithelial_initial <- make_epithelial_counts(all_human, "Initial Dataset")
epithelial_hq <- all_human %>%
  semi_join(hq_epithelial_mags, by = "MAG") %>%
  make_epithelial_counts("HQ MAGs (90%, 5%)")
epithelial_final <- all_human %>%
  filter(taxonomy_status == "final_reannotation") %>%
  make_epithelial_counts("After Selection")

# =========================
# 3. Environmental habitats
# =========================
environment_initial <- tribble(
  ~Stage, ~Classification, ~Source, ~Value,
  "Initial Dataset", "Environmental", "Aquatic", 19300,
  "Initial Dataset", "Environmental", "Other environmental", 5941,
  "Initial Dataset", "Environmental", "Terrestrial", 3410
)

environment_hq <- hq_data %>%
  filter(
    str_starts(check_Source, "Aquatic") |
      str_starts(check_Source, "other") |
      str_starts(check_Source, "Terrestrial")
  ) %>%
  mutate(
    Source = case_when(
      str_starts(check_Source, "Aquatic") ~ "Aquatic",
      str_starts(check_Source, "Terrestrial") ~ "Terrestrial",
      TRUE ~ "Other environmental"
    )
  ) %>%
  count(Source, name = "Value") %>%
  mutate(Stage = "HQ MAGs (90%, 5%)", Classification = "Environmental") %>%
  select(Stage, Classification, Source, Value)

environment_final <- final_metadata %>%
  filter(from == "gem") %>%
  mutate(Source = recode(source, other = "Other environmental")) %>%
  count(Source, name = "Value") %>%
  mutate(Stage = "After Selection", Classification = "Environmental") %>%
  select(Stage, Classification, Source, Value)

# =========================
# 4. Reconciliation and plotting data
# =========================
stage_levels <- c("Initial Dataset", "HQ MAGs (90%, 5%)", "After Selection")
source_levels <- c(
  "Dry skin",
  "Sebaceous skin",
  "Moist skin",
  "Toenail",
  "Nasal cavity",
  "Unknown",
  "Terrestrial",
  "Aquatic",
  "Other environmental"
)

source_colors <- c(
  "Dry skin" = "#B2182B",
  "Sebaceous skin" = "#D6604D",
  "Moist skin" = "#F4A582",
  "Toenail" = "#FDDBC7",
  "Nasal cavity" = "#C51B7D",
  "Unknown" = "#BDBDBD",
  "Terrestrial" = "#1B7837",
  "Aquatic" = "#2166AC",
  "Other environmental" = "#67A9CF"
)

selection_data <- bind_rows(
  epithelial_initial,
  epithelial_hq,
  epithelial_final,
  environment_initial,
  environment_hq,
  environment_final
) %>%
  mutate(
    Stage = factor(Stage, levels = stage_levels),
    Classification = factor(Classification, levels = c("Epithelial-associated", "Environmental")),
    Source = factor(Source, levels = source_levels)
  )

expected_totals <- tribble(
  ~Stage, ~Classification, ~Expected,
  "Initial Dataset", "Epithelial-associated", 22748,
  "Initial Dataset", "Environmental", 28651,
  "HQ MAGs (90%, 5%)", "Epithelial-associated", 7262,
  "HQ MAGs (90%, 5%)", "Environmental", 11869,
  "After Selection", "Epithelial-associated", 3811,
  "After Selection", "Environmental", 3437
)

observed_totals <- selection_data %>%
  mutate(Stage = as.character(Stage), Classification = as.character(Classification)) %>%
  group_by(Stage, Classification) %>%
  summarise(Observed = sum(Value), .groups = "drop")

reconciliation <- expected_totals %>%
  left_join(observed_totals, by = c("Stage", "Classification"))

stopifnot(
  nrow(all_human) == 22748,
  n_distinct(all_human$MAG) == 22748,
  nrow(hq_data) == 19131,
  nrow(final_metadata) == 7248,
  nrow(hq_epithelial_mags) == 7262,
  sum(all_human$taxonomy_status == "final_reannotation") == 3811,
  all(reconciliation$Expected == reconciliation$Observed),
  sum(selection_data$Value[selection_data$Stage == "Initial Dataset"]) == 51399,
  sum(selection_data$Value[selection_data$Stage == "HQ MAGs (90%, 5%)"]) == 19131,
  sum(selection_data$Value[selection_data$Stage == "After Selection"]) == 7248
)

plot_data <- selection_data %>%
  complete(Stage, Classification, Source, fill = list(Value = 0)) %>%
  arrange(Stage, Classification, Source) %>%
  group_by(Stage, Classification) %>%
  mutate(
    Stage_Index = as.numeric(Stage),
    X = Stage_Index + if_else(Classification == "Epithelial-associated", -0.20, 0.20),
    Ymax = cumsum(Value),
    Ymin = Ymax - Value,
    Ymid = (Ymin + Ymax) / 2,
    Label = if_else(Value >= 300, scales::comma(Value), ""),
    Label_Color = if_else(Source %in% c("Moist skin", "Toenail", "Unknown", "Other environmental"), "#303030", "white"),
    Xmin = X - 0.17,
    Xmax = X + 0.17
  ) %>%
  ungroup() %>%
  filter(Value > 0)

totals <- plot_data %>%
  group_by(Stage, Classification, X) %>%
  summarise(Total = sum(Value), .groups = "drop")

class_labels <- totals %>%
  mutate(Label = if_else(Classification == "Epithelial-associated", "EPITHELIAL\nASSOCIATED", "ENVIRONMENTAL\n "))

subtitle_text <- "Overall retention: 14.1% | Epithelial-associated: 16.8% | Environmental: 12.0%"

# =========================
# 5. Main panel and grouped right-side legend
# =========================
p_main <- ggplot(plot_data) +
  geom_rect(
    aes(xmin = Xmin, xmax = Xmax, ymin = Ymin, ymax = Ymax, fill = Source),
    color = "white",
    linewidth = 0.28,
    alpha = 0.98
  ) +
  geom_text(
    aes(x = X, y = Ymid, label = Label, color = Label_Color),
    size = 2.75,
    fontface = "bold",
    show.legend = FALSE
  ) +
  geom_label(
    data = totals,
    aes(x = X, y = Total, label = scales::comma(Total)),
    inherit.aes = FALSE,
    vjust = -0.42,
    size = 3.05,
    fontface = "bold",
    label.padding = unit(0.14, "lines"),
    linewidth = 0.15,
    fill = "white",
    color = "grey15"
  ) +
  geom_text(
    data = class_labels,
    aes(x = X, y = 0, label = Label),
    inherit.aes = FALSE,
    vjust = 1.45,
    size = 2.15,
    fontface = "bold",
    color = "grey30"
  ) +
  scale_fill_manual(values = source_colors, breaks = source_levels, drop = FALSE) +
  scale_color_identity() +
  scale_x_continuous(
    breaks = seq_along(stage_levels),
    labels = stage_levels,
    expand = expansion(mult = c(0.08, 0.08))
  ) +
  scale_y_continuous(
    labels = scales::comma,
    expand = expansion(mult = c(0.07, 0.14))
  ) +
  labs(
    title = "Genome retention across curation stages",
    subtitle = subtitle_text,
    x = "Curation stage",
    y = "Number of MAGs"
  ) +
  coord_cartesian(clip = "off") +
  theme_classic(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 13.5, hjust = 0, color = "grey15", margin = margin(b = 5)),
    plot.subtitle = element_text(size = 8.7, color = "grey35", hjust = 0, margin = margin(b = 14)),
    axis.line = element_line(color = "grey45", linewidth = 0.35),
    axis.ticks = element_line(color = "grey45", linewidth = 0.3),
    axis.text.x = element_text(face = "bold", size = 8.9, color = "grey25"),
    axis.text.y = element_text(size = 8.5, color = "grey30"),
    axis.title = element_text(face = "bold", size = 9.5, color = "grey20"),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.28),
    legend.position = "none",
    plot.margin = margin(16, 4, 50, 16)
  )

legend_items <- tribble(
  ~kind, ~label, ~source, ~y,
  "header", "Epithelial-associated sites", NA_character_, 9.5,
  "item", "Dry skin", "Dry skin", 8.65,
  "item", "Sebaceous skin", "Sebaceous skin", 8.05,
  "item", "Moist skin", "Moist skin", 7.45,
  "item", "Toenail", "Toenail", 6.85,
  "item", "Nasal cavity", "Nasal cavity", 6.25,
  "item", "Unknown", "Unknown", 5.65,
  "header", "Environmental habitats", NA_character_, 4.45,
  "item", "Terrestrial", "Terrestrial", 3.60,
  "item", "Aquatic", "Aquatic", 3.00,
  "item", "Other environmental", "Other environmental", 2.40
) %>%
  mutate(
    source = factor(source, levels = source_levels),
    color = source_colors[as.character(source)]
  )

p_legend <- ggplot() +
  geom_text(
    data = filter(legend_items, kind == "header"),
    aes(x = 0, y = y, label = label),
    hjust = 0,
    size = 3.15,
    fontface = "bold",
    color = "grey18"
  ) +
  geom_tile(
    data = filter(legend_items, kind == "item"),
    aes(x = 0.08, y = y, fill = source),
    width = 0.18,
    height = 0.24,
    color = "white",
    linewidth = 0.25
  ) +
  geom_text(
    data = filter(legend_items, kind == "item"),
    aes(x = 0.22, y = y, label = label),
    hjust = 0,
    size = 2.9,
    color = "grey25"
  ) +
  scale_fill_manual(values = source_colors, limits = source_levels, drop = FALSE) +
  coord_cartesian(xlim = c(0, 2.25), ylim = c(1.85, 9.85), clip = "off") +
  theme_void() +
  theme(legend.position = "none", plot.margin = margin(20, 8, 50, 0))

footnote_text <- paste0(
  "Epithelial-associated MAGs are grouped by body site or skin microenvironment; nasal cavity is shown as a body site.\n",
  "Other environmental habitats include built environment, wastewater, solid waste, bioremediation, industrial production, and air."
)

# =========================
# 6. Save output
# =========================
dest_dir <- "outputs/Figure1/"
dir.create(dest_dir, recursive = TRUE, showWarnings = FALSE)
write_tsv(plot_data, paste0(dest_dir, "Figure1B_Selection_Efficiency_Modern.source_data.tsv"))

draw_figure <- function() {
  grid.newpage()
  print(p_main, vp = viewport(x = 0.00, y = 0.10, width = 0.74, height = 0.88, just = c("left", "bottom")))
  print(p_legend, vp = viewport(x = 0.745, y = 0.18, width = 0.245, height = 0.76, just = c("left", "bottom")))
  grid.text(
    footnote_text,
    x = unit(0.035, "npc"),
    y = unit(0.03, "npc"),
    just = c("left", "bottom"),
    gp = gpar(fontsize = 7.0, col = "grey35")
  )
}

pdf(paste0(dest_dir, "Figure1B_Selection_Efficiency_Modern.pdf"), width = 11.0, height = 7.8, useDingbats = FALSE)
draw_figure()
dev.off()

png(paste0(dest_dir, "Figure1B_Selection_Efficiency_Modern.png"), width = 3300, height = 2340, res = 300, type = "cairo")
draw_figure()
dev.off()

print(reconciliation)
print(plot_data %>% arrange(Stage, Classification, Source))
print("Single-panel Figure 1B saved as PDF and PNG.")
