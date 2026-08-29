library(ggplot2)
library(dplyr)
library(tidyr)

# Paths
input_path <- "data/figure3/Grand_Butterfly_AllRanks_Data.tsv"
output_dir <- "outputs/Figure3"

if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# Load and basic cleaning
df <- read.table(input_path, sep = "\t", header = TRUE, check.names = FALSE)

# Remove the CF22/CF25 panel from Figure S2B.
excluded_cfs <- c("CF22", "CF25")
df <- df %>%
  filter(!CF_A %in% excluded_cfs, !CF_B %in% excluded_cfs)

df$Rank <- factor(df$Rank, levels = c("Phylum", "Class", "Order", "Family", "Genus"))

# 1. Rename and reorder 'Other'
df <- df %>%
  mutate(Color_Group = ifelse(Color_Group == "Exclusion/Other", "Other", Color_Group))

# Get unique phyla that are not "Other"
colored_phyla <- sort(unique(df$Color_Group[df$Color_Group != "Other"]))

# 2. Setup user-defined color palette
user_colors <- c(
  "#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C",
  "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A", "#FFFF99", "#B15928"
)

# Map colors to phyla (cycling if phyla count > 12, though unlikely here)
phylum_mapping <- setNames(rep(user_colors, length.out = length(colored_phyla)), colored_phyla)

# Add "Other" as Grey at the end
final_palette <- c(phylum_mapping, "Other" = "#D3D3D3")

# Ensure Color_Group is a factor with 'Other' at the end for the legend
df$Color_Group <- factor(df$Color_Group, levels = c(colored_phyla, "Other"))

# 3. Plot: The Grand Butterfly Map (v8.2 - User Style)
p <- ggplot(df, aes(group = Taxon)) +
  # 1. Connectivity Segment (The Wings)
  geom_segment(aes(x = Copy_A, y = Prev_A, xend = Copy_B, yend = Prev_B, color = Color_Group),
    alpha = 0.4, size = 0.4
  ) +
  # 2. Endpoint A (Left)
  geom_point(aes(x = Copy_A, y = Prev_A, color = Color_Group, size = N_Total), alpha = 0.7) +
  # 3. Endpoint B (Right)
  geom_point(aes(x = Copy_B, y = Prev_B, color = Color_Group, size = N_Total), alpha = 0.7) +
  # 4. Facet Grid 5x5
  facet_grid(Rank ~ Pattern, scales = "free_x") +
  # 5. Mirror X Axis (Absolute labels)
  scale_x_continuous(labels = abs) +
  scale_y_continuous(limits = c(0, 115), breaks = c(0, 50, 100)) +
  # Apply User Palette
  scale_color_manual(values = final_palette, name = "Phylum (Colored if both subfamilies present)") +
  scale_size_continuous(range = c(0.8, 4), guide = "none") +
  # 6. Central Axis
  geom_vline(xintercept = 0, color = "black", linetype = "solid", size = 0.6) +
  # Styling
  theme_bw() +
  labs(
    title = "Figure 3D: Hierarchical Implementation of ecological Candidate Factor Subfamilies",
    subtitle = "Segments connect A/B implementations of the same lineage. Colored = Co-occurrence | Grey = Lineage-specific locking.",
    x = "Mean Copy Number",
    y = "Occurrence Frequency (%)"
  ) +
  theme(
    strip.background = element_rect(fill = "gray15"),
    strip.text = element_text(color = "white", face = "bold", size = 9),
    panel.grid.minor = element_blank(),
    axis.text = element_text(size = 7),
    legend.position = "bottom",
    panel.spacing = unit(0.5, "lines")
  )

# 4. Correct Subfamily Labels (Using edge-positioning to avoid overlap across different X scales)
df_labs <- df %>%
  group_by(Pattern) %>%
  summarise(CF_A = unique(CF_A), CF_B = unique(CF_B), .groups = "drop") %>%
  pivot_longer(cols = c(CF_A, CF_B), names_to = "Side", values_to = "Label") %>%
  mutate(
    x_pos = ifelse(Side == "CF_A", -Inf, Inf), 
    y_pos = 112,
    hjust_val = ifelse(Side == "CF_A", -0.1, 1.1) # Pad away from edges
  )

p <- p + geom_text(
  data = df_labs, aes(x = x_pos, y = y_pos, label = sub("^CF", "eCF", Label), hjust = hjust_val),
  size = 3.5, fontface = "bold", color = "black", inherit.aes = FALSE
)

# Save
ggsave(file.path(output_dir, "Figure3D_Grand_Butterfly_Map_v8.2.pdf"), p, width = 18, height = 15)
ggsave(file.path(output_dir, "Figure3D_Grand_Butterfly_Map_v8.2.png"), p, width = 18, height = 15, dpi = 300)

print(paste("Refined Grand Map generated successfully in:", output_dir))
