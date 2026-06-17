# =========================
# 0. Load Packages
# =========================
library(tidyverse)
library(patchwork)

# =========================
# 1. Read Data
# =========================
data <- read_tsv("E:/data/filter/data/03_selected_genomes_warning_priority.tsv")

# =========================
# 2. Habitat Definition
# =========================
data <- data %>%
  mutate(
    Habitat = case_when(
      from == "gem" ~ "Environmental",
      from == "human" & source == "nose" ~ "Nasal",
      from == "human" ~ "Skin",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(Habitat))

# =========================
# 3. Taxonomy Truncation Function (Core)
# =========================
extract_lineage_level <- function(taxonomy, level) {
  
  levels <- c("d", "p", "c", "o", "f", "g", "s")
  idx <- match(level, levels)
  
  # Split taxonomy
  parts <- str_split(taxonomy, ";")[[1]]
  
  # If level is not enough, pad with NA
  if (length(parts) < idx) {
    return(NA)
  }
  
  # Truncate to corresponding level
  lineage <- paste(parts[1:idx], collapse = ";")
  
  # ===== p / c suffix merge =====
  if (level %in% c("p", "c")) {
    lineage <- str_replace(
      lineage,
      paste0("(", level, "__[^;]+)_[A-Z]$"),
      "\\1"
    )
  }
  
  return(lineage)
}

# =========================
# 4. Extract all levels
# =========================
processed_data <- data %>%
  mutate(
    Phylum  = map_chr(classification, ~extract_lineage_level(.x, "p")),
    Class   = map_chr(classification, ~extract_lineage_level(.x, "c")),
    Order   = map_chr(classification, ~extract_lineage_level(.x, "o")),
    Family  = map_chr(classification, ~extract_lineage_level(.x, "f")),
    Genus   = map_chr(classification, ~extract_lineage_level(.x, "g")),
    Species = map_chr(classification, ~extract_lineage_level(.x, "s"))
  )

# =========================
# 5. Diversity Statistics (unique count per level)
# =========================
diversity_summary <- processed_data %>%
  group_by(Habitat) %>%
  summarise(
    Phylum  = n_distinct(Phylum, na.rm = TRUE),
    Class   = n_distinct(Class, na.rm = TRUE),
    Order   = n_distinct(Order, na.rm = TRUE),
    Family  = n_distinct(Family, na.rm = TRUE),
    Genus   = n_distinct(Genus, na.rm = TRUE),
    Species = n_distinct(Species, na.rm = TRUE)
  ) %>%
  pivot_longer(
    cols = -Habitat,
    names_to = "Taxonomic_Level",
    values_to = "Count"
  )

# =========================
# 6. Sort Levels (ensure correct order)
# =========================
diversity_summary$Taxonomic_Level <- factor(
  diversity_summary$Taxonomic_Level,
  levels = c("Phylum", "Class", "Order", "Family", "Genus", "Species")
)
# =========================
# 7. Plotting (facet display + independent y-axis + right legend + text labels)
# =========================
p <- ggplot(diversity_summary, 
            aes(x = Habitat, 
                y = Count, 
                fill = Habitat)) +
  
  # Bar chart
  geom_bar(stat = "identity", 
           width = 0.7,
           color = "white", 
           linewidth = 0.2) +
  
  # Core: Add numerical labels
  geom_text(aes(label = scales::comma(Count)), # Use comma-separated number format
            vjust = -0.5,                      # Labels above bars
            size = 3.5,                        # Font size
            fontface = "bold") +               # Bold labels
  
  # Facets
  facet_wrap(~Taxonomic_Level, 
             scales = "free_y", 
             nrow = 2) + 
  
  # Axis expansion: increase mult upper limit to ensure labels are not clipped
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.2)), 
    labels = scales::comma
  ) +
  
  scale_fill_manual(values = c(
    "Environmental" = "#3288BD",
    "Nasal" = "#FEE08B",
    "Skin" = "#D53E4F"
  )) +
  
  labs(
    x = NULL,
    y = "Number of taxa",
    fill = "Habitat"
  ) +
  
  theme_bw(base_size = 12) + 
  theme(
    strip.background = element_rect(fill = "grey95", color = "grey80"),
    strip.text = element_text(face = "bold", size = 10),
    axis.text.x = element_blank(),     
    axis.ticks.x = element_blank(),
    panel.grid.major.x = element_blank(),
    legend.position = "right", 
    legend.title = element_text(face = "bold"),
    plot.margin = margin(10, 10, 10, 10)
  )

print(p)
# =========================
# 8. Save
# =========================
ggsave("Figure1C_Taxonomic_Diversity.pdf",
       p, width = 8, height = 6, dpi = 300)

ggsave("Figure1C_Taxonomic_Diversity.png",
       p, width = 8, height = 6, dpi = 300)
