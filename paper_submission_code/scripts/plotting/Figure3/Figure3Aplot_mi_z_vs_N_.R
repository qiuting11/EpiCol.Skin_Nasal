library(ggplot2)
library(dplyr)
library(ggrepel)

# 1. Load configuration
source("R/gephe_R/scripts/Figure3_Final/Figure3_color_config.R")

# 2. Read data
data_path <- "gephe_output_11/analysis/output_test_1/top0.2_I3.0_modules.txt"
df <- read.delim(data_path, stringsAsFactors = FALSE)

# 3. Data mapping
df <- df %>%
  mutate(
    theme = functional_theme_map[module],
    module = factor(module, levels = module_order)
  )

# 4. Plotting: Theme fill + Module border
p <- ggplot(df, aes(x = N, y = mi_z, fill = theme, color = module)) +
  geom_point(size = 4, shape = 21, stroke = 1.5, alpha = 0.9) +
  scale_fill_manual(values = functional_theme_colors, name = "Functional Theme") +
  scale_color_manual(values = full_module_colors, name = "Module") +
  scale_x_log10(breaks = c(1, 10, 100, 1000, 3000)) + 
  geom_text_repel(data = subset(df, N >= 50 | mi_z >= 900), aes(label = cf_id), color = "black", size = 3) +
  theme_bw(base_size = 14) +
  guides(fill = guide_legend(override.aes = list(shape = 21, stroke = 0.5, color = "grey30")),
         color = guide_legend(override.aes = list(shape = 21, fill = "white"))) +
  labs(x = "Quantity (N, log10 scale)", y = "MI Z-score", title = "V4: Theme (Fill) & Module (Border)")

ggsave("gephe_output_11/analysis/output_test_1/MI_Zscore_vs_N_v4.pdf", p, width = 11, height = 6.5)
