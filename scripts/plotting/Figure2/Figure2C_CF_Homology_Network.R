library(tidyverse)
library(igraph)
library(RColorBrewer)

base_path <- "data/external/figure2_homology_network/"
required_files <- file.path(base_path, c("pog_relationships.tsv", "pog_info.tsv"))
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0) {
  stop(
    "Missing external Figure2C input(s): ", paste(missing_files, collapse = ", "),
    "\nSee docs/data_availability.md for the expected external data layout."
  )
}

pog_rels <- read_delim(file.path(base_path, "pog_relationships.tsv"), delim = "\t", show_col_types = FALSE)
pog_info <- read_delim(file.path(base_path, "pog_info.tsv"), delim = "\t", show_col_types = FALSE)
modules <- read_delim("data/figure2/top0.2_I3.0_modules.txt", delim = "\t", show_col_types = FALSE)

pog_cf_mapping <- modules %>%
  mutate(pog_num = as.numeric(str_replace(original_pog, "pog", ""))) %>%
  select(pog_num, cf_id) %>%
  distinct()

dat_graph <- pog_rels %>%
  mutate(link = link / (nrow(pog_info)^2)) %>%
  mutate(weight = log(link + 1)) %>%
  select(from = pog1, to = pog2, weight) %>%
  mutate(across(c(from, to), as.character))

graph <- graph_from_data_frame(dat_graph, directed = FALSE, vertices = as.character(pog_info$pog)) %>%
  igraph::simplify(edge.attr.comb = "sum")

V(graph)$label <- pog_cf_mapping$cf_id[match(as.numeric(V(graph)$name), pog_cf_mapping$pog_num)]
V(graph)$label <- sub("^CF", "eCF", V(graph)$label)
V(graph)$label <- ifelse(is.na(V(graph)$label), V(graph)$name, V(graph)$label)

V(graph)$color <- ifelse(
  as.numeric(V(graph)$name) %in% pog_cf_mapping$pog_num,
  adjustcolor("SkyBlue2", alpha.f = .5),
  adjustcolor("gray", alpha.f = .3)
)

my_draw_logic <- function() {
  plot(graph,
    layout = layout_with_dh(graph),
    vertex.size = 12,
    vertex.frame.color = "white",
    vertex.label.color = "black",
    vertex.label.cex = 0.9,
    vertex.label.dist = 1.5,
    vertex.label.degree = -pi / 2,
    edge.width = 2,
    edge.color = adjustcolor("gray40", alpha.f = 0.5),
    margin = -0.05
  )
}

dir.create("outputs/Figure2", recursive = TRUE, showWarnings = FALSE)
png(file.path("outputs/Figure2", "Figure2D_CF_Homology_Network.png"), width = 2000, height = 2000, res = 300)
my_draw_logic()
dev.off()

pdf(file.path("outputs/Figure2", "Figure2D_CF_Homology_Network.pdf"), width = 8, height = 8)
my_draw_logic()
dev.off()
