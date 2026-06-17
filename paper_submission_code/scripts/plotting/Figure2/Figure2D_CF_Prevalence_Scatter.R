# Load necessary packages
library(ggplot2)
library(dplyr)
library(ggrepel)

# Read data
setwd("E:/data/R/gephe_R/results/figures/archives/figure2_legacy_results")
top_data <- read.csv("2C_fast/top0.2_statistics_with_modules.csv")
random_data <- read.csv("2c_output_1/2C_random_statistics_with_preference.csv")

# Use the provided color scheme
module_colors <- c(
  "#72190E", "#f94144", "#F7557F", "#7C417F", "#A05992",
  "#f3722c", "#f8961e", "#f9c74f", "#F7F056",
  "#C2B923", "#A9D88C", "#90be6d", "#43aa8b",
  "#277da1", "#1BA3C6", "#B9F2F0", "#4F7CBA", "#5289C7"
)

# Sort by module number order (instead of alphabetical order)
# Extract module numbers and sort
module_order <- unique(top_data$module)
module_nums <- as.numeric(gsub("module", "", module_order))
sorted_modules <- module_order[order(module_nums)]
names(module_colors) <- sorted_modules

# Reorder factor levels in the data frame to ensure the legend follows numerical order
top_data$module <- factor(top_data$module, levels = sorted_modules)

# Calculate the number of POGs for each module to be used in the legend
module_counts <- top_data %>%
  group_by(module) %>%
  summarise(n_pogs = n(), .groups = "drop") %>%
  arrange(factor(module, levels = sorted_modules))

# Create legend labels: module name + number of POGs
module_labels <- paste0(module_counts$module, " (n=", module_counts$n_pogs, ")")
names(module_labels) <- module_counts$module

# Create main plot
p <- ggplot() +
  # Background layer: random dataset (gray, low transparency)
  geom_point(
    data = random_data,
    aes(x = gem_freq, y = human_freq),
    color = "gray85", alpha = 0.25, size = 1.5, shape = 16
  ) +

  # Foreground layer: Top 0.2% dataset colored by module - using solid points
  geom_point(
    data = top_data,
    aes(x = gem_freq, y = human_freq, color = module),
    size = 4, alpha = 0.95, shape = 16
  ) + # shape=16 is a solid circle

  # Add diagonal line: represents equal frequencies
  geom_abline(
    intercept = 0, slope = 1, linetype = "dashed",
    color = "black", alpha = 0.6, linewidth = 0.8
  ) +

  # Annotation layer: label each point with cf_id
  geom_text_repel(
    data = top_data,
    aes(x = gem_freq, y = human_freq, label = cf_id),
    size = 3.5, box.padding = 0.35, point.padding = 0.25,
    min.segment.length = 0.25, max.overlaps = 30,
    color = "black", fontface = "bold",
    segment.color = "gray40", segment.size = 0.3
  ) +

  # Set coordinate axes - 1:1 ratio
  scale_x_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, 0.2),
    labels = scales::percent_format(accuracy = 1),
    expand = expansion(mult = c(0.02, 0.05))
  ) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, 0.2),
    labels = scales::percent_format(accuracy = 1),
    expand = expansion(mult = c(0.02, 0.05))
  ) +

  # Set colors and legend - using solid points
  scale_color_manual(
    name = "Functional Modules",
    values = module_colors,
    labels = module_labels,
    guide = guide_legend(
      override.aes = list(size = 3, alpha = 1, shape = 16) # Ensure legend points are solid
    )
  ) +

  # Labels and titles
  labs(
    x = "Environmental Frequency (n = 3,437 MAGs)",
    y = "Host-associated Frequency (n = 3,811 MAGs)",
    title = "Prevalence of Candidate Colonization Factors",
  ) +

  # Ensure 1:1 ratio
  coord_fixed(ratio = 1) +

  # Theme settings
  theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(
      face = "bold",
      size = 16,
      hjust = 0.5,
      margin = margin(b = 5)
    ),
    # Axis labels
    axis.title = element_text(face = "bold", size = 12),
    axis.text = element_text(size = 10, color = "black"),


    # Legend - ensure numerical order
    legend.position = "right",
    legend.title = element_text(
      face = "bold",
      size = 11,
      margin = margin(b = 5)
    ),
    legend.text = element_text(size = 7),
    legend.key.size = unit(0.6, "cm"),
    legend.spacing.y = unit(0.1, "cm"),

    # Grid and background
    panel.grid.major = element_line(
      color = "gray88",
      linewidth = 0.4,
      linetype = "solid"
    ),
    panel.grid.minor = element_line(
      color = "gray92",
      linewidth = 0.2,
      linetype = "solid"
    ),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    panel.border = element_rect(fill = NA, color = "gray70", linewidth = 0.5),

    # Margins
    plot.margin = margin(20, 30, 20, 20)
  )

# Display plot
print(p)

# Save high-resolution images
ggsave("E:/data/R/gephe_R/results/figures/figure2/FigureS2C_CF_Prevalence_Scatter.pdf",
  plot = p,
  width = 10,
  height = 8,
  dpi = 300
)

ggsave("E:/data/R/gephe_R/results/figures/figure2/FigureS2C_CF_Prevalence_Scatter.png",
  plot = p,
  width = 10,
  height = 8,
  dpi = 300
)

cat("Plot saved successfully!\n")
cat("1. FigureS2C_CF_Prevalence_Scatter.pdf - 1:1 ratio corrected version\n")
cat("2. FigureS2C_CF_Prevalence_Scatter.png - 1:1 ratio corrected version\n")
