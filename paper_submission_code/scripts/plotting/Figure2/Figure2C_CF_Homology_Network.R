library(tidyverse)
library(igraph)
library(RColorBrewer)

# 1. Unified data reading
base_path <- "E:/data/gephe_output_11/analysis/output_test_1/"
pog_rels <- read_delim(paste0(base_path, "pog_relationships.tsv"), delim = "\t", show_col_types = FALSE)
pog_info <- read_delim(paste0(base_path, "pog_info.tsv"), delim = "\t", show_col_types = FALSE)
modules  <- read_delim(paste0(base_path, "top0.2_I3.0_modules.txt"), delim = "\t", show_col_types = FALSE)

# 2. Process mapping and edge data (merged into one pipeline)
pog_cf_mapping <- modules %>%
  mutate(pog_num = as.numeric(str_replace(original_pog, "pog", ""))) %>%
  select(pog_num, cf_id) %>%
  distinct()

dat_graph <- pog_rels %>%
  mutate(link = link / (nrow(pog_info)^2)) %>%
  mutate(weight = log(link + 1)) %>% # Perform log transformation directly here
  select(from = pog1, to = pog2, weight) %>%
  mutate(across(c(from, to), as.character))

# 3. Construct and simplify the network (in one step)
graph <- graph_from_data_frame(dat_graph, directed = FALSE, vertices = as.character(pog_info$pog)) %>%
  igraph::simplify(edge.attr.comb = "sum") # Merge duplicate edges and sum weights

# 4. Prepare visual attributes
V(graph)$label  <- pog_cf_mapping$cf_id[match(as.numeric(V(graph)$name), pog_cf_mapping$pog_num)]
V(graph)$label  <- ifelse(is.na(V(graph)$label), V(graph)$name, V(graph)$label)

V(graph)$color  <- ifelse(as.numeric(V(graph)$name) %in% pog_cf_mapping$pog_num, 
                          adjustcolor("SkyBlue2", alpha.f = .5), 
                          adjustcolor("gray", alpha.f = .3))

# Define a plotting function
my_draw_logic <- function() {
  plot(graph,
       layout = layout_with_dh(graph),
       vertex.size = 12,
       vertex.frame.color = "white",
       vertex.label.color = 'black',
       vertex.label.cex = 0.9,
       vertex.label.dist = 1.5,
       vertex.label.degree = -pi/2,
       edge.width = 2, # Thicken the lines slightly
       edge.color = adjustcolor("gray40", alpha.f = 0.5),
       margin = -0.05)
}

# Save PNG
png("Figure2D_CF_Homology_Network.png", width=2000, height=2000, res=300); my_draw_logic(); dev.off()

# Save PDF
pdf("Figure2D_CF_Homology_Network.pdf", width=8, height=8); my_draw_logic(); dev.off()
