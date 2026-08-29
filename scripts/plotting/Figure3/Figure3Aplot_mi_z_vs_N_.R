library(ggplot2)
library(dplyr)
library(ggrepel)

source("scripts/plotting/Figure3/Figure3_color_config.R")

df <- read.delim("data/figure3/top0.2_I3.0_modules.txt", stringsAsFactors = FALSE)

df <- df %>%
  mutate(
    theme = functional_theme_map[module],
    module = factor(module, levels = module_order)
  )

p <- ggplot(df, aes(x = N, y = mi_z, fill = theme, color = module)) +
  geom_point(size = 4, shape = 21, stroke = 1.5, alpha = 0.9) +
  scale_fill_manual(values = functional_theme_colors, name = "Functional theme") +
  scale_color_manual(values = full_module_colors, labels = sub("^module", "eCM", names(full_module_colors)), name = "eCM") +
  scale_x_log10(breaks = c(1, 10, 100, 1000, 3000)) +
  geom_text_repel(data = subset(df, N >= 50 | mi_z >= 900), aes(label = sub("^CF", "eCF", cf_id)), color = "black", size = 3) +
  theme_bw(base_size = 14) +
  guides(
    fill = guide_legend(override.aes = list(shape = 21, stroke = 0.5, color = "grey30")),
    color = guide_legend(override.aes = list(shape = 21, fill = "white"))
  ) +
  labs(x = "Quantity (N, log10 scale)", y = "MI Z-score", title = "MI Z-score vs. eCF prevalence")

dir.create("outputs/Figure3", recursive = TRUE, showWarnings = FALSE)
ggsave(file.path("outputs/Figure3", "MI_Zscore_vs_N_v4.pdf"), p, width = 11, height = 6.5)
