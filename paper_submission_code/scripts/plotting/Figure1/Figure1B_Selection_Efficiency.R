library(tidyverse)

# =========================
# 1. Data Definition
# =========================
selection_data <- tribble(
  ~Stage, ~Classification, ~Source, ~Value,
  "Initial Dataset", "Human", "HSMG", 15995,
  "Initial Dataset", "Human", "UHSG", 5779,
  "Initial Dataset", "Human", "Nasal", 974,
  "Initial Dataset", "Environmental", "Aquatic", 19300,
  "Initial Dataset", "Environmental", "Other", 5941,
  "Initial Dataset", "Environmental", "Terrestrial", 3410,
  "HQ MAGs(90%,5%)", "Human", "HSMG", 4471,
  "HQ MAGs(90%,5%)", "Human", "UHSG", 2112,
  "HQ MAGs(90%,5%)", "Human", "Nasal", 679,
  "HQ MAGs(90%,5%)", "Environmental", "Aquatic", 7088,
  "HQ MAGs(90%,5%)", "Environmental", "Other", 3326,
  "HQ MAGs(90%,5%)", "Environmental", "Terrestrial", 1455,
  "After Selection", "Human", "HSMG", 2553,
  "After Selection", "Human", "UHSG", 928,
  "After Selection", "Human", "Nasal", 330,
  "After Selection", "Environmental", "Aquatic", 1958,
  "After Selection", "Environmental", "Other", 1095,
  "After Selection", "Environmental", "Terrestrial", 384
) %>%
  mutate(
    Stage = factor(Stage, levels = c("Initial Dataset", "HQ MAGs(90%,5%)", "After Selection")),
    Classification = factor(Classification, levels = c("Human", "Environmental")),
    Source = factor(Source, levels = c("HSMG", "UHSG", "Nasal", "Aquatic", "Other", "Terrestrial"))
  )

# Calculate stats for the subtitle
initial_total <- sum(selection_data$Value[selection_data$Stage == "Initial Dataset"])
final_total <- sum(selection_data$Value[selection_data$Stage == "After Selection"])
overall_ratio <- (final_total / initial_total) * 100

human_initial <- sum(selection_data$Value[selection_data$Stage == "Initial Dataset" & selection_data$Classification == "Human"])
human_final <- sum(selection_data$Value[selection_data$Stage == "After Selection" & selection_data$Classification == "Human"])
human_ratio <- (human_final / human_initial) * 100

env_initial <- sum(selection_data$Value[selection_data$Stage == "Initial Dataset" & selection_data$Classification == "Environmental"])
env_final <- sum(selection_data$Value[selection_data$Stage == "After Selection" & selection_data$Classification == "Environmental"])
env_ratio <- (env_final / env_initial) * 100

subtitle_text <- sprintf(
  "Overall Selection Ratio: %.1f%% | Human Sources: %.1f%% | Environmental Sources: %.1f%%",
  overall_ratio, human_ratio, env_ratio
)

# Colors matching Python script but polished for R
# We add dummy entries for legend headers
source_colors <- c(
  "HUMAN" = "white",
  "HSMG" = "#D62728",
  "UHSG" = "#FF7F0E",
  "Nasal" = "#FFBB78",
  "GEM" = "white",
  "Aquatic" = "#1F77B4",
  "Other" = "#2CA02C",
  "Terrestrial" = "#17BECF"
)

# Define the legend order (Top-down)
legend_breaks <- c("HUMAN", "Nasal", "UHSG", "HSMG", "GEM", "Terrestrial", "Other", "Aquatic")

# Calculate totals for labels
totals <- selection_data %>%
  group_by(Stage, Classification) %>%
  summarise(Total = sum(Value), .groups = "drop")

# =========================
# 2. Plotting
# =========================
# Ensure Source order matches stacking (from bottom to top)
# And include dummy headers in levels
selection_data <- selection_data %>%
  mutate(Source = factor(Source, levels = c("HSMG", "UHSG", "Nasal", "Aquatic", "Other", "Terrestrial", "HUMAN", "GEM")))

p <- ggplot(selection_data, aes(x = Classification, y = Value, fill = Source)) +
  # Stacked bars (default position_stack is bottom-to-top)
  geom_col(position = "stack", width = 0.7, color = "white", linewidth = 0.2) +
  # Faceting by Stage
  facet_grid(~Stage, switch = "x") +
  # Value labels inside bars (only if > 300)
  geom_text(
    aes(label = ifelse(Value > 300, scales::comma(Value), ""), 
        color = ifelse(Value > 1500, "white", "black")),
    position = position_stack(vjust = 0.5),
    size = 3, fontface = "bold", show.legend = FALSE
  ) +
  # Total labels on top
  geom_label(
    data = totals,
    aes(x = Classification, y = Total, label = scales::comma(Total)),
    vjust = -0.5, size = 3.5, fontface = "bold",
    label.padding = unit(0.2, "lines"),
    fill = "white", alpha = 0.8,
    inherit.aes = FALSE
  ) +
  # Color and scale customizations
  scale_fill_manual(
    values = source_colors,
    breaks = legend_breaks,
    drop = FALSE, # Ensure dummy levels stay in legend
    guide = guide_legend(
      override.aes = list(
        # Hide the box for headers by setting alpha or color to white/transparent
        fill = c("white", source_colors[c("Nasal", "UHSG", "HSMG")], "white", source_colors[c("Terrestrial", "Other", "Aquatic")]),
        color = c("white", "white", "white", "white", "white", "white", "white", "white")
      )
    )
  ) +
  scale_color_identity() +
  scale_y_continuous(labels = scales::comma, expand = expansion(mult = c(0, 0.15))) +
  # Labels and Title
  labs(
    title = "Genome Selection Efficiency",
    subtitle = subtitle_text,
    x = NULL,
    y = "Number of MAGs",
    fill = NULL
  ) +
  # Theme customization
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5, margin = margin(b = 10)),
    plot.subtitle = element_text(size = 10, color = "grey20", hjust = 0.5, margin = margin(b = 15)),
    axis.title.y = element_text(face = "bold"),
    axis.text.x = element_text(face = "bold", color = "grey30"),
    strip.text = element_text(face = "bold", size = 10),
    strip.placement = "outside",
    panel.grid.major.x = element_blank(),
    legend.position = "right",
    legend.text = element_text(size = 9),
    plot.margin = margin(15, 15, 15, 15)
  )

# =========================
# 3. Save Output
# =========================
ggsave("E:/data/R/gephe_R/scripts/Figure1/Figure1B_Selection_Efficiency.pdf",
  p,
  width = 9, height = 7, dpi = 300
)
ggsave("E:/data/R/gephe_R/scripts/Figure1/Figure1B_Selection_Efficiency.png",
  p,
  width = 8.5, height = 7, dpi = 300
)

print("Visualization saved to E:/data/R/gephe_R/scripts/Figure1/Figure1B_Selection_Efficiency.png")
