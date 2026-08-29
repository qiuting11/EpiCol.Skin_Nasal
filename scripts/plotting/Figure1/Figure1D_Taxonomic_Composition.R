library(tidyverse)
library(RColorBrewer)
library(patchwork)

# =========================
# =========================
data <- read_tsv("data/figure1/03_selected_genomes_warning_priority.tsv")

# =========================
# =========================
data <- data %>%
  mutate(
    Habitat = case_when(
      from == "gem" ~ "Environmental",
      from == "human" ~ "Epithelial-associated",
      TRUE ~ NA_character_
    )
  )

# =========================
# =========================
process_taxonomy <- function(taxonomy_string, level) {
  pattern <- paste0(";", level, "__([^;]+)")
  matches <- str_extract(taxonomy_string, pattern)
  
  if (is.na(matches)) return(NA)
  
  taxon <- str_replace(matches, paste0(";", level, "__"), "")
  
  if (level %in% c("p", "c")) {
    taxon <- str_replace(taxon, "_[A-Z]$", "")
  }
  
  return(taxon)
}

# =========================
# =========================
processed_data <- data %>%
  mutate(
    Phylum = map_chr(classification, ~process_taxonomy(.x, "p")),
    Class  = map_chr(classification, ~process_taxonomy(.x, "c"))
  ) %>%
  filter(!is.na(Habitat), !is.na(Phylum), !is.na(Class)) %>%
  select(Habitat, Phylum, Class)

# =========================
# =========================
prepare_plot_data <- function(df, tax_level, top_n = 10) {
  
  stats <- df %>%
    count(Habitat, !!sym(tax_level)) %>%
    group_by(Habitat) %>%
    mutate(Percentage = n / sum(n) * 100) %>%
    ungroup()
  
  top_taxa <- stats %>%
    group_by(Habitat) %>%
    slice_max(n, n = top_n) %>%
    pull(!!sym(tax_level)) %>%
    unique()
  
  stats <- stats %>%
    mutate(
      group = ifelse(!!sym(tax_level) %in% top_taxa,
                     !!sym(tax_level), "Other")
    )
  
  plot_data <- stats %>%
    group_by(Habitat, group) %>%
    summarise(n = sum(n), .groups = "drop") %>%
    group_by(Habitat) %>%
    mutate(Percentage = n / sum(n) * 100) %>%
    ungroup()
  
  global_order <- plot_data %>%
    group_by(group) %>%
    summarise(total = sum(n), .groups = "drop") %>%
    arrange(desc(total)) %>%
    pull(group)
  
  plot_data <- plot_data %>%
    mutate(
      group = factor(group, levels = c(setdiff(global_order, "Other"), "Other"))
    ) %>%
    group_by(Habitat) %>%
    arrange(group) %>%
    mutate(
      ypos = cumsum(Percentage) - 0.5 * Percentage
    ) %>%
    ungroup()
  
  return(plot_data)
}

# =========================
# =========================
phylum_plot_data <- prepare_plot_data(processed_data, "Phylum", 10)
class_plot_data  <- prepare_plot_data(processed_data, "Class", 10)

# =========================
# =========================
custom_colors <- c(
  "#A6CEE3", "#579CC7", "#3688AD", "#8BC395", "#89CB6C", 
  "#40A635", "#919D5F", "#F99392", "#EB494A", "#E83C2D",
  "#F79C5D", "#FDA746", "#FE8205", "#E39970", "#BFA5CF",
  "#8861AC", "#917099", "#E7E099", "#DEB969", "#B15928"
)

other_color <- "#808080"

get_colors <- function(groups) {
  
  unique_groups <- levels(groups)
  n <- length(unique_groups) - 1
  
  colors <- rep(custom_colors, length.out = n)
  colors <- c(colors, other_color)
  
  names(colors) <- unique_groups
  
  return(colors)
}

# =========================
# =========================
plot_bar <- function(plot_data, legend_name) {
  
  colors <- get_colors(plot_data$group)
  
  ggplot(plot_data, aes(x = Habitat, y = Percentage, fill = group)) +
    geom_bar(stat = "identity", width = 0.7, color = "white", linewidth = 0.2) +
    scale_fill_manual(values = colors, name = legend_name) +
    labs(
      x = NULL,
      y = "Percentage (%)"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      axis.text.x = element_text(face = "bold"),
      panel.grid.major.x = element_blank(),
      legend.position = "right",
      legend.key.height = unit(0.4, "cm"),
      legend.text = element_text(size = 7),
      legend.title = element_text(size = 9, face = "bold")
    )
}

# =========================
# =========================
phylum_plot <- plot_bar(phylum_plot_data, "Phylum") +
  guides(fill = guide_legend(ncol = 1))

class_plot  <- plot_bar(class_plot_data, "Class") +
  guides(fill = guide_legend(ncol = 1))

combined_plot <- phylum_plot + class_plot +
  plot_layout(ncol = 2, widths = c(1, 1)) &
  theme(
    plot.margin = margin(3, 5, 3, 5)
  )

print(combined_plot)

# =========================
# =========================
ggsave(file.path("outputs/Figure1", "Figure1D_Taxonomic_Composition.pdf"),
       combined_plot, width = 10, height = 8, dpi = 300)

ggsave(file.path("outputs/Figure1", "Figure1D_Taxonomic_Composition.png"),
       combined_plot, width = 10, height = 8, dpi = 300)
