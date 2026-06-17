library(tidyverse)
library(patchwork)

# =========================
# 1. Raw Data (51k)
# =========================
gem_raw  <- read_tsv("E:/data/checkm2/gem_all_quality_reports.tsv")
human_raw <- read_tsv("E:/data/checkm2/skin_nose_all_quality_reports.tsv")

raw_data <- bind_rows(gem_raw, human_raw) %>%
  mutate(
    Dataset = "Raw",
    Group = case_when(
      str_detect(form, "nose|skin") ~ "Human",
      TRUE ~ "Environmental"
    )
  )

# =========================
# 2. Filtered (19k)
# =========================
filtered_data <- read_tsv("E:/data/filter/data/02_merged_data_with_sources_categories.tsv") %>%
  mutate(
    Dataset = "Filtered",
    Group = case_when(
      from == "human" ~ "Human",
      from == "gem" ~ "Environmental"
    )
  )

# =========================
# 3. Final (7k)
# =========================
final_data <- read_tsv("E:/data/filter/data/03_selected_genomes_warning_priority.tsv") %>%
  mutate(
    Dataset = "Final",
    Group = case_when(
      from == "human" ~ "Human",
      from == "gem" ~ "Environmental"
    )
  )

# =========================
# 4. Merge
# =========================
combined_data <- bind_rows(raw_data, filtered_data, final_data) %>%
  filter(!is.na(Group)) %>%
  mutate(
    Dataset = factor(Dataset, levels = c("Raw", "Filtered", "Final")),
    Group = factor(Group, levels = c("Human", "Environmental"))
  )

# =========================
# 5. Plotting Function
# =========================
create_plot <- function(df, y_var, y_label, log_scale = FALSE, ylim = NULL) {
  
  p <- ggplot(df, aes(x = Group, y = .data[[y_var]], fill = Dataset)) +
    geom_violin(position = position_dodge(0.8), width = 0.7,
                alpha = 0.6, color = NA, trim = TRUE) +
    geom_boxplot(position = position_dodge(0.8), width = 0.15,
                 outlier.shape = NA, alpha = 0.9) +
    scale_fill_manual(values = c(
      "Raw" = "#579CC7",
      "Filtered" ="#89CB6C",
      "Final" = "#E83C2D"
    )) +
    labs(x = NULL, y = y_label) +
    theme_minimal(base_size = 11) +
    theme(
      axis.text.x = element_text(face = "bold"),
      panel.grid.major.x = element_blank(),
      legend.position = "bottom",
      legend.title = element_blank()
    )
  
  if (log_scale) {
    p <- p + scale_y_log10(expand = c(0, 0))
  } else {
    p <- p + scale_y_continuous(expand = expansion(mult = c(0, 0.02)))
  }
  
  if (!is.null(ylim)) {
    p <- p + coord_cartesian(ylim = ylim)
  }
  
  return(p)
}

# =========================
# 6. Three Metrics
# =========================

# Completeness
p1 <- create_plot(combined_data, "Completeness", "Completeness (%)", ylim = c(25, 100)) +
  geom_hline(yintercept = 90, linetype = "dashed", color = "grey40")

p2 <- create_plot(combined_data, "Contamination", "Contamination (%)", ylim = c(0, 10)) +
  geom_hline(yintercept = 5, linetype = "dashed", color = "grey40")

p3 <- create_plot(combined_data, "Contig_N50", "Contig N50 (bp)", log_scale = TRUE)

# =========================
# 7. Layout (1x3)
# =========================
final_plot <- p1 + p2 + p3 +
  plot_layout(ncol = 3, guides = "collect") &
  theme(
    legend.position = "bottom"
  )

print(final_plot)

# =========================
# 8. Save
# =========================
ggsave("FigureS1_Quality_ThreeStages.pdf",
       final_plot, width = 11, height = 4, dpi = 300)

ggsave("FigureS1_Quality_ThreeStages.png",
       final_plot, width = 11, height = 4, dpi = 300)
