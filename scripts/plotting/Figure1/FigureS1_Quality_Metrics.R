library(tidyverse)
library(patchwork)

# =========================
# 1. Source files
# =========================
gem_raw_file <- "data/figure1/gem_all_quality_reports.tsv"
epithelial_raw_file <- "data/figure1/skin_nose_all_quality_reports.tsv"
final_file <- "data/figure1/03_selected_genomes_warning_priority.tsv"

# =========================
# 2. Raw and final datasets
# =========================
gem_raw <- read_tsv(gem_raw_file, show_col_types = FALSE) %>%
  mutate(Dataset = "Raw", Group = "Environmental")

epithelial_raw <- read_tsv(epithelial_raw_file, show_col_types = FALSE) %>%
  mutate(Dataset = "Raw", Group = "Epithelial-associated")

final_data <- read_tsv(final_file, show_col_types = FALSE) %>%
  mutate(
    Dataset = "Final",
    Group = case_when(
      from == "human" ~ "Epithelial-associated",
      from == "gem" ~ "Environmental",
      TRUE ~ NA_character_
    )
  )

combined_data <- bind_rows(gem_raw, epithelial_raw, final_data) %>%
  filter(!is.na(Group)) %>%
  mutate(
    Dataset = factor(Dataset, levels = c("Raw", "Final")),
    Group = factor(Group, levels = c("Epithelial-associated", "Environmental"))
  )

# Keep a compact source table for auditability.
source_data <- combined_data %>%
  select(Dataset, Group, Completeness, Contamination, Contig_N50)

# =========================
# 3. Reconciliation
# =========================
expected_counts <- tribble(
  ~Dataset, ~Group, ~Expected,
  "Raw", "Epithelial-associated", 22748,
  "Raw", "Environmental", 28651,
  "Final", "Epithelial-associated", 3811,
  "Final", "Environmental", 3437
)

observed_counts <- combined_data %>%
  mutate(Dataset = as.character(Dataset), Group = as.character(Group)) %>%
  count(Dataset, Group, name = "Observed")

reconciliation <- expected_counts %>%
  left_join(observed_counts, by = c("Dataset", "Group"))

stopifnot(
  nrow(gem_raw) == 28651,
  nrow(epithelial_raw) == 22748,
  nrow(final_data) == 7248,
  all(reconciliation$Expected == reconciliation$Observed),
  all(c("Completeness", "Contamination", "Contig_N50") %in% colnames(combined_data))
)

# =========================
# 4. Plot function
# =========================
stage_colors <- c(
  "Raw" = "#4C72B0",
  "Final" = "#DD8452"
)

dodge <- position_dodge(width = 0.78)

create_metric_plot <- function(data, y_var, y_label, title, log_scale = FALSE,
                               y_limits = NULL, threshold = NULL,
                               threshold_label = NULL, threshold_y = NULL) {
  p <- ggplot(data, aes(x = Group, y = .data[[y_var]], fill = Dataset)) +
    geom_violin(
      alpha = 0.58,
      width = 0.76,
      color = "white",
      linewidth = 0.18,
      trim = TRUE,
      scale = "width",
      position = dodge
    ) +
    geom_boxplot(
      width = 0.14,
      alpha = 0.88,
      outlier.shape = NA,
      linewidth = 0.28,
      position = dodge
    ) +
    stat_summary(
      fun = median,
      geom = "point",
      aes(group = Dataset),
      position = dodge,
      shape = 21,
      size = 1.9,
      stroke = 0.35,
      fill = "white",
      color = "black"
    ) +
    stat_summary(
      fun = mean,
      geom = "point",
      aes(group = Dataset),
      position = dodge,
      shape = 23,
      size = 1.9,
      stroke = 0.35,
      fill = "#D62728",
      color = "black"
    ) +
    scale_fill_manual(values = stage_colors, name = "Dataset stage") +
    labs(title = title, x = NULL, y = y_label) +
    theme_minimal(base_size = 10) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 10.5, margin = margin(b = 5)),
      axis.title = element_text(face = "bold", size = 9.5),
      axis.text.x = element_text(face = "bold", size = 8.6, color = "grey25"),
      axis.text.y = element_text(size = 8.2, color = "grey30"),
      legend.position = "bottom",
      legend.title = element_text(face = "bold", size = 9),
      legend.text = element_text(size = 8.5),
      panel.grid.major.x = element_blank(),
      panel.grid.major.y = element_line(color = "grey90", linewidth = 0.28),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "grey82", fill = NA, linewidth = 0.45),
      plot.margin = margin(5, 6, 5, 6)
    )

  if (log_scale) {
    p <- p +
      scale_y_log10(labels = scales::label_comma()) +
      annotation_logticks(sides = "l", linewidth = 0.25)
  } else {
    p <- p + scale_y_continuous(labels = scales::label_number())
  }

  if (!is.null(y_limits)) {
    p <- p + coord_cartesian(ylim = y_limits)
  }

  if (!is.null(threshold)) {
    p <- p + geom_hline(yintercept = threshold, linetype = "dashed", color = "grey35", linewidth = 0.35)
  }

  if (!is.null(threshold_label)) {
    label_y <- ifelse(is.null(threshold_y), threshold, threshold_y)
    p <- p + annotate("text", x = 0.62, y = label_y, label = threshold_label,
                      color = "grey25", size = 2.6, hjust = 0)
  }

  p
}

# =========================
# 5. Metrics
# =========================
p1 <- create_metric_plot(
  source_data,
  "Completeness",
  "Completeness (%)",
  "Genome completeness",
  y_limits = c(20, 100),
  threshold = 90,
  threshold_label = "90% threshold",
  threshold_y = 92
)

p2 <- create_metric_plot(
  source_data,
  "Contamination",
  "Contamination (%)",
  "Genome contamination",
  y_limits = c(0, 12),
  threshold = 5,
  threshold_label = "5% threshold",
  threshold_y = 5.8
)

p3 <- create_metric_plot(
  source_data,
  "Contig_N50",
  "Contig N50 (bp)",
  "Assembly contiguity (N50)",
  log_scale = TRUE
)

final_plot <- p1 + p2 + p3 +
  plot_layout(ncol = 3, guides = "collect") &
  theme(legend.position = "bottom")

final_plot <- final_plot +
  plot_annotation(
    title = "Quality metrics before and after genome curation",
    subtitle = "Raw and final MAG datasets grouped as epithelial-associated or environmental",
    caption = paste(
      "Raw: epithelial-associated n=22,748, environmental n=28,651;",
      "Final: epithelial-associated n=3,811, environmental n=3,437.",
      "Violin plots show distributions; boxes show quartiles; white circles show medians; red diamonds show means."
    ),
    theme = theme(
      plot.title = element_text(size = 13, face = "bold", hjust = 0.5, margin = margin(b = 6)),
      plot.subtitle = element_text(size = 9.5, hjust = 0.5, color = "grey35", margin = margin(b = 12)),
      plot.caption = element_text(size = 7.5, hjust = 0, color = "grey45", lineheight = 1.15, margin = margin(t = 8))
    )
  )

print(final_plot)

# =========================
# 6. Save output
# =========================
dest_dir <- "outputs/Figure1/"
dir.create(dest_dir, recursive = TRUE, showWarnings = FALSE)
write_tsv(source_data, paste0(dest_dir, "FigureS1_Quality_Metrics.source_data.tsv"))
write_tsv(reconciliation, paste0(dest_dir, "FigureS1_Quality_Metrics.reconciliation.tsv"))

ggsave(
  paste0(dest_dir, "FigureS1_Quality_Metrics.pdf"),
  final_plot,
  width = 11,
  height = 4.6,
  dpi = 300,
  device = cairo_pdf
)

ggsave(
  paste0(dest_dir, "FigureS1_Quality_Metrics.png"),
  final_plot,
  width = 11,
  height = 4.6,
  dpi = 300,
  bg = "white"
)

print(reconciliation)
print("Figure S1 quality metrics saved as PDF and PNG.")
