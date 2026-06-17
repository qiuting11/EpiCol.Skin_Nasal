library(tidyverse)

# =========================
# 1. Data Preparation
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

# Calculate stage-to-stage retention ratios for annotation
retention <- selection_data %>%
  group_by(Stage) %>%
  summarise(Total = sum(Value)) %>%
  mutate(
    Prev_Total = lag(Total),
    Ratio = Total / Prev_Total * 100,
    Label = sprintf("Retention: %.1f%%", Ratio)
  ) %>%
  filter(!is.na(Ratio))

# Calculate classification totals
totals <- selection_data %>%
  group_by(Stage, Classification) %>%
  summarise(Total = sum(Value), .groups = "drop")

# Colors: Specified Spectral-like Palette
# HSMG/UHSG (Skin) are reddish-orange
# Nasal is yellow
# Aquatic/Terrestrial are greens
# Other is grey
modern_colors <- c(
  "HSMG"        = "#D73027", # Deep Red (Skin sub)
  "UHSG"        = "#FDAE61", # Orange (Skin sub/Requested)
  "Nasal"       = "#FEE08B", # Yellow (Nasal)
  "Aquatic"     = "#ABDDA4", # Light Green
  "Terrestrial" = "#66C2A5", # Green
  "Other"       = "grey80"   # Light Grey
)

# =========================
# 2. Modern Plotting
# =========================
p <- ggplot(selection_data, aes(x = Classification, y = Value, fill = Source)) +
  # Use a slight transparency for a softer look
  geom_col(width = 0.7, alpha = 0.9, color = "white", linewidth = 0.2) +
  # Use free_y to ensure the "After Selection" bars are not too small
  # ncol = 2 will force 3 stages into 2 rows
  facet_wrap(~Stage, scales = "free_y", strip.position = "top", ncol = 2) +
  
  # Value labels inside bars (centered)
  geom_text(
    aes(label = ifelse(Value > 300, scales::comma(Value), ""), 
        color = ifelse(Value > 5000 | Source == "HSMG" | Source == "Aquatic", "white", "black")),
    position = position_stack(vjust = 0.5),
    size = 3.2, fontface = "bold", show.legend = FALSE
  ) +
  
  # Total labels on top of each bar
  geom_text(
    data = totals,
    aes(x = Classification, y = Total, label = scales::comma(Total)),
    inherit.aes = FALSE,
    vjust = -0.5, size = 3.8, fontface = "bold", color = "grey10"
  ) +
  
  scale_fill_manual(values = modern_colors) +
  scale_color_identity() +
  scale_y_continuous(labels = scales::comma, expand = expansion(mult = c(0, 0.18))) +
  
  labs(
    title = "GENOME SELECTION EFFICIENCY",
    subtitle = "Relative composition and scaling across selection stages (Y-axes are scaled independently)",
    x = NULL, y = "Number of Genomes",
    fill = "Source Environment"
  ) +
  
  theme_minimal(base_size = 13) +
  theme(
    # Titles
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5, color = "#2C3E50"),
    plot.subtitle = element_text(size = 10, color = "grey40", hjust = 0.5, margin = margin(b = 20)),
    
    # Strip (Stage Labels) - Clean & Bold
    strip.background = element_rect(fill = "grey98", color = NA),
    strip.text = element_text(face = "bold", size = 11, color = "grey20", margin = margin(t=5, b=5)),
    
    # Grid & Axis
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(color = "grey92", linetype = "dashed"),
    panel.spacing = unit(2, "lines"), # Increase spacing between rows
    
    axis.line.x = element_line(color = "grey70"),
    axis.text.x = element_text(face = "bold", size = 10, color = "grey30"),
    axis.title.y = element_text(face = "bold", margin = margin(r = 10)),
    
    # Legend
    legend.position = "bottom",
    legend.title = element_text(size = 10, face = "bold"),
    legend.text = element_text(size = 9),
    
    plot.margin = margin(20, 20, 20, 20)
  )

# =========================
# 3. Save Output
# =========================
dest_dir <- "E:/data/R/gephe_R/scripts/Figure1/"
# Adjust height for 2-row layout
ggsave(paste0(dest_dir, "Figure1B_Selection_Efficiency_Modern.pdf"), p, width = 9, height = 11)
ggsave(paste0(dest_dir, "Figure1B_Selection_Efficiency_Modern.png"), p, width = 9, height = 11, dpi = 300)

print("Modern version saved as Figure1B_Selection_Efficiency_Modern.png")
