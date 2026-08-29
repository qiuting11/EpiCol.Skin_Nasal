library(gggenes)
library(ggplot2)
library(dplyr)

input_file <- "data/figure5/metNIQ_absolute_synteny_data_fixed.tsv"
output_dir <- "outputs/Figure5"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

data <- read.delim(input_file, header = TRUE, sep = "\t")

data$gene_label <- case_when(
  data$gene == "metQ" ~ "metQ (Binding)",
  data$gene == "metI" ~ "metI (Permease)",
  data$gene == "metN" ~ "metN (ATPase)",
  TRUE ~ data$gene
)
data$gene_label <- factor(data$gene_label, levels = c("metQ (Binding)", "metI (Permease)", "metN (ATPase)"))

met_palette <- c("#E69F00", "#56B4E9", "#009E73")

p <- ggplot(data, aes(xmin = start, xmax = end, y = molecule, fill = gene_label, forward = (strand == "forward"))) +
  geom_gene_arrow(arrowhead_height = unit(4, "mm"), arrowhead_width = unit(1, "mm"), size = 0.5) +
  facet_wrap(~molecule, scales = "free", ncol = 1, strip.position = "left") +
  scale_fill_manual(values = met_palette) +
  theme_genes() +
  theme(
    panel.spacing.y = unit(0.2, "lines"),
    strip.text.y.left = element_text(angle = 0, face = "italic", size = 9, hjust = 1),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_text(size = 8),
    legend.position = "bottom",
    legend.title = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    strip.background = element_blank(),
    panel.border = element_blank(),
    axis.line.x = element_line(color = "black")
  ) +
  labs(
    title = "metNIQ Operon: Real Genomic Coordinates",
    x = "Genomic Position on Contig (bp)"
  )

ggsave(file.path(output_dir, "Figure5_metNIQ_real_coordinates.png"), plot = p, width = 10, height = 12, dpi = 300)
ggsave(file.path(output_dir, "Figure5_metNIQ_real_coordinates.pdf"), plot = p, width = 10, height = 12)

message("Figure5 metNIQ synteny plot saved to outputs/Figure5")
