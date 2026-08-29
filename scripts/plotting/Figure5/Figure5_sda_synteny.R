library(gggenes)
library(ggplot2)
library(dplyr)

input_file <- "data/figure5/sda_synteny_clean_core_only_scaled.tsv"
output_dir <- "outputs/Figure5"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

data <- read.delim(input_file, header = TRUE, sep = "\t")

# Filter out blank rows so the x-axis tightens around genes.
plot_data <- data %>% filter(type != "blank")

mol_order <- c(
  "Peptococcus niger (human2946)",
  "Enterococcus faecalis (human1001)",
  "Staphylococcus hominis (human1)",
  "Cutibacterium modestum (human102)",
  "Corynebacterium kefirresidentii (human0)",
  "Pseudomonas vlassakiae (human246)"
)
plot_data$molecule <- factor(plot_data$molecule, levels = mol_order)

sda_palette <- c(
  "sdaA (Fused)" = "#acccf6",
  "sdaAA (Alpha)" = "#db5682",
  "sdaAB (Beta)" = "#f1131e"
)

p <- ggplot(plot_data, aes(xmin = start, xmax = end, y = molecule, fill = type, forward = (strand == "forward"))) +
  geom_gene_arrow(arrowhead_height = unit(4, "mm"), arrowhead_width = unit(1, "mm"), size = 0.5, color = "black") +
  facet_wrap(~molecule, scales = "free", ncol = 1, strip.position = "left") +
  scale_fill_manual(values = sda_palette) +
  scale_x_continuous(expand = expansion(mult = c(0.1, 0.1))) +
  theme_genes() +
  theme(
    panel.spacing.y = unit(0.2, "lines"),
    strip.text.y.left = element_text(angle = 0, face = "italic", size = 9, hjust = 1),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_text(size = 8, color = "black"),
    axis.line.x = element_line(color = "black", linewidth = 0.5),
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(size = 9),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    strip.background = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank()
  ) +
  labs(
    title = "Structural Comparison of Serine Dehydratase Architectures",
    x = "Genomic Position (bp)"
  )

ggsave(file.path(output_dir, "Figure5_sda_synteny.png"), plot = p, width = 10, height = 7, dpi = 300)
ggsave(file.path(output_dir, "Figure5_sda_synteny.pdf"), plot = p, width = 10, height = 7)

message("Figure5 SDA synteny plot saved to outputs/Figure5")
