library(ggplot2)
library(dplyr)
library(ggrepel)
# Input files are resolved relative to repository root; no setwd needed.
top_data <- read.csv("data/figure2/top0.2_statistics_with_modules.csv")
random_data <- read.csv("data/figure2/2C_random_statistics_with_preference.csv")

module_colors <- c(
  "#72190E", "#f94144", "#F7557F", "#7C417F", "#A05992",
  "#f3722c", "#f8961e", "#f9c74f", "#F7F056",
  "#C2B923", "#A9D88C", "#90be6d", "#43aa8b",
  "#277da1", "#1BA3C6", "#B9F2F0", "#4F7CBA", "#5289C7"
)

module_order <- unique(top_data$module)
module_nums <- as.numeric(gsub("module", "", module_order))
sorted_modules <- module_order[order(module_nums)]
names(module_colors) <- sorted_modules

top_data$module <- factor(top_data$module, levels = sorted_modules)

module_counts <- top_data %>%
  group_by(module) %>%
  summarise(n_pogs = n(), .groups = "drop") %>%
  arrange(factor(module, levels = sorted_modules))

module_labels <- paste0(sub("^module", "eCM", module_counts$module), " (n=", module_counts$n_pogs, ")")
names(module_labels) <- module_counts$module

p <- ggplot() +
  geom_point(
    data = random_data,
    aes(x = gem_freq, y = human_freq),
    color = "gray85", alpha = 0.25, size = 1.5, shape = 16
  ) +

  geom_point(
    data = top_data,
    aes(x = gem_freq, y = human_freq, color = module),
    size = 4, alpha = 0.95, shape = 16
  ) +

  geom_abline(
    intercept = 0, slope = 1, linetype = "dashed",
    color = "black", alpha = 0.6, linewidth = 0.8
  ) +

  geom_text_repel(
    data = top_data,
    aes(x = gem_freq, y = human_freq, label = sub("^CF", "eCF", cf_id)),
    size = 3.5, box.padding = 0.35, point.padding = 0.25,
    min.segment.length = 0.25, max.overlaps = 30,
    color = "black", fontface = "bold",
    segment.color = "gray40", segment.size = 0.3
  ) +

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

  scale_color_manual(
    name = "ecological Candidate Modules",
    values = module_colors,
    labels = module_labels,
    guide = guide_legend(
      override.aes = list(size = 3, alpha = 1, shape = 16)
    )
  ) +

  labs(
    x = "Environmental Frequency (n = 3,437 MAGs)",
    y = "Epithelial-associated Frequency (n = 3,811 MAGs)",
    title = "Prevalence of ecological Candidate Factors",
  ) +

  coord_fixed(ratio = 1) +

  theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(
      face = "bold",
      size = 16,
      hjust = 0.5,
      margin = margin(b = 5)
    ),
    axis.title = element_text(face = "bold", size = 12),
    axis.text = element_text(size = 10, color = "black"),


    legend.position = "right",
    legend.title = element_text(
      face = "bold",
      size = 11,
      margin = margin(b = 5)
    ),
    legend.text = element_text(size = 7),
    legend.key.size = unit(0.6, "cm"),
    legend.spacing.y = unit(0.1, "cm"),

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

    plot.margin = margin(20, 30, 20, 20)
  )

print(p)

ggsave(file.path("outputs/Figure2", "FigureS2C_CF_Prevalence_Scatter.pdf"),
  plot = p,
  width = 10,
  height = 8,
  dpi = 300
)

ggsave(file.path("outputs/Figure2", "FigureS2C_CF_Prevalence_Scatter.png"),
  plot = p,
  width = 10,
  height = 8,
  dpi = 300
)

cat("Figure2D prevalence scatter completed.\n")
cat("Outputs: outputs/Figure2/FigureS2C_CF_Prevalence_Scatter.pdf/png\n")
