# ==============================================================================
# ==============================================================================
if (!requireNamespace("gtools", quietly = TRUE)) install.packages("gtools")
library(dplyr)
library(tidyr)
library(ComplexHeatmap)
library(circlize)
library(gtools)
library(grid)
ht_opt$message <- FALSE

# ==============================================================================
# ==============================================================================

modules_df <- read.delim("data/figure2/top0.2_I3.0_modules.txt", sep = "\t")
metadata <- read.delim("data/figure1/metadata_for_analysis.tsv", sep = "\t")

pog_matrix <- as.matrix(read.table("data/figure2/pog_phylogenetic_profiles.tsv",
  sep = "\t", header = TRUE, row.names = 1, check.names = FALSE
))
storage.mode(pog_matrix) <- "numeric"

cf_pog_mapping <- modules_df %>%
  mutate(pog_num = as.numeric(gsub("pog", "", original_pog))) %>%
  select(cf_id, pog_num, original_pog) %>%
  distinct()

cf_genome_mat <- pog_matrix[paste0("POG_", cf_pog_mapping$pog_num), , drop = FALSE]
rownames(cf_genome_mat) <- cf_pog_mapping$cf_id

cf_sim_mat_raw <- read.table("data/figure2/top0.2_I3.0_cocluster.txt",
  header = TRUE, row.names = 1, check.names = FALSE
)
cf_sim_mat <- as.matrix(cf_sim_mat_raw)
new_names <- cf_pog_mapping$cf_id[match(rownames(cf_sim_mat), cf_pog_mapping$original_pog)]
rownames(cf_sim_mat) <- colnames(cf_sim_mat) <- new_names

common_cfs <- intersect(rownames(cf_genome_mat), rownames(cf_sim_mat))
cf_genome_mat <- cf_genome_mat[common_cfs, ]
cf_sim_mat <- cf_sim_mat[common_cfs, common_cfs]

# ==============================================================================
# ==============================================================================

all_modules <- unique(modules_df$module)
sorted_module_levels <- gtools::mixedsort(all_modules)
sorted_module_labels <- sub("^module", "eCM", sorted_module_levels)

annot_module_df <- modules_df %>%
  select(cf_id, module) %>%
  distinct()
annot_module_df <- annot_module_df[match(common_cfs, annot_module_df$cf_id), ]
annot_module_df$module <- factor(sub("^module", "eCM", annot_module_df$module), levels = sorted_module_labels)

my_colors <- c(
  "#72190E", "#f94144", "#F7557F", "#7C417F", "#A05992", "#f3722c", "#f8961e",
  "#f9c74f", "#F7F056", "#C2B923", "#A9D88C", "#90be6d", "#43aa8b", "#277da1",
  "#1BA3C6", "#B9F2F0", "#4F7CBA", "#5289C7"
)
module_colors <- structure(my_colors[1:length(sorted_module_labels)], names = sorted_module_labels)

metadata_to_use <- data.frame(taxon_oid = colnames(cf_genome_mat)) %>%
  left_join(metadata, by = "taxon_oid")

phylum_counts <- table(metadata_to_use$Phylum)
top20_phyla <- names(sort(phylum_counts, decreasing = TRUE))[1:20]

metadata_to_use <- metadata_to_use %>%
  mutate(
    Habitat = factor(case_when(
      source %in% c("HSMG", "UHSG") ~ "Skin",
      source == "nose" ~ "Nose",
      source %in% c("Aquatic", "Terrestrial", "other") ~ source,
      TRUE ~ "Unknown"
    ), levels = c("Skin", "Nose", "Aquatic", "Terrestrial", "other", "Unknown")),
    From = factor(case_when(
      from == "human" ~ "Epithelial-associated",
      from == "gem" ~ "Environmental",
      TRUE ~ "Unknown"
    ), levels = c("Epithelial-associated", "Environmental", "Unknown")),
    Phylum_Plot = ifelse(Phylum %in% top20_phyla, Phylum, "Other"),
    Phylum_Plot = factor(Phylum_Plot, levels = c(top20_phyla, "Other"))
  )

phylum_palette <- colorRampPalette(c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33", "#A65628", "#F781BF"))(21)
names(phylum_palette) <- levels(metadata_to_use$Phylum_Plot)
phylum_palette["Other"] <- "grey80"

# ==============================================================================
# ==============================================================================

annot_col_genome <- HeatmapAnnotation(
  Phylum = metadata_to_use$Phylum_Plot,
  Habitat = metadata_to_use$Habitat,
  Origin = metadata_to_use$From,
  col = list(
    Phylum = phylum_palette,
    Habitat = c("Skin" = "#FDAE61", "Nose" = "#FEE08B", "Aquatic" = "#ABDDA4", "Terrestrial" = "#66C2A5", "other" = "grey80", "Unknown" = "grey40"),
    Origin = c("Environmental" = "#3288BD", "Epithelial-associated" = "#D53E4F", "Unknown" = "grey40")
  ),
  annotation_legend_param = list(
    Phylum = list(ncol = 1, title_position = "topleft"),
    Habitat = list(ncol = 1, title_position = "topleft"),
    Origin = list(ncol = 1, title_position = "topleft")
  )
)

ha_row_module <- rowAnnotation(
  Module = annot_module_df$module,
  col = list(Module = module_colors),
  annotation_legend_param = list(
    Module = list(ncol = 1, title_position = "topleft", at = sorted_module_labels)
  ),
  show_annotation_name = FALSE
)

# ==============================================================================
# ==============================================================================

col_sim <- colorRamp2(c(0, 0.5, 1), c("white", "#FFE5B4", "#E41A1C"))
col_presence <- colorRamp2(c(0, 1), c("white", "grey"))

ht1 <- Heatmap(
  cf_sim_mat,
  name = "Similarity",
  col = col_sim,
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  clustering_distance_rows = function(x) as.dist(1 - x),
  clustering_distance_columns = function(x) as.dist(1 - x),
  show_row_names = FALSE,
  show_column_names = FALSE,
  width = unit(20, "cm"),
  heatmap_legend_param = list(ncol = 1, title_position = "topleft")
)

ht2 <- Heatmap(
  cf_genome_mat,
  name = "Presence",
  col = col_presence,
  cluster_rows = FALSE,
  cluster_columns = TRUE,
  top_annotation = annot_col_genome,
  right_annotation = ha_row_module,
  show_row_names = TRUE,
  row_labels = sub("^CF", "eCF", rownames(cf_genome_mat)),
  row_names_side = "right",
  show_column_names = FALSE,
  column_title = paste0(ncol(cf_genome_mat), " Genomes"),
  use_raster = TRUE,
  raster_quality = 5,
  heatmap_legend_param = list(ncol = 1, title_position = "topleft")
)

# ==============================================================================
# ==============================================================================

pdf(file.path("outputs/Figure2", "FigureS2E_CF_Coinheritance_PhyloProfile_ByPhylum.pdf"), width = 20, height = 8.7)
draw(ht1 + ht2,
  main_heatmap = "Similarity", heatmap_legend_side = "right",
  annotation_legend_side = "right", legend_grouping = "original"
)
dev.off()

png(file.path("outputs/Figure2", "FigureS2E_CF_Coinheritance_PhyloProfile_ByPhylum.png"), width = 20, height = 8, units = "in", res = 300)
draw(ht1 + ht2,
  main_heatmap = "Similarity", heatmap_legend_side = "right",
  annotation_legend_side = "right", legend_grouping = "original"
)
dev.off()

cat("FigureS2E eCF/eCM phylo-profile by phylum saved.\n")
