# Load necessary libraries  
library(dplyr)
library(ggalluvial)
library(ggplot2)
# Read data
mcl_1 <- read.csv('data/figure2/mcl_1.4.csv')
mcl_2 <- read.csv('data/figure2/mcl_2.0.csv')
mcl_3 <- read.csv('data/figure2/mcl_3.0.csv')

# Rename columns to avoid conflicts
colnames(mcl_1) <- c("protein", "pog_I1.4")
colnames(mcl_2) <- c("protein", "pog_I2.0") 
colnames(mcl_3) <- c("protein", "pog_I3.0")

# Merge data
dat <- mcl_1 %>%
  left_join(mcl_2, by = 'protein') %>%
  left_join(mcl_3, by = 'protein') %>%
  group_by(pog_I1.4, pog_I2.0, pog_I3.0) %>%
  dplyr::summarize(size = n(), .groups = 'drop')

my_colors <- c("#A6CEE3", "#76AFD2", "#4690C1", "#287EB1", "#5CA3A2", 
               "#90C793", "#A1D67D", "#74C05C", "#47AA3B", "#599E41", 
               "#A09C67", "#E79A8E", "#F47878", "#EC4B4C", "#E31E20", 
               "#EB4F36", "#F48954", "#FDBC6B", "#FDA644", "#FE8F1C", 
               "#F98314", "#E79660", "#D4A8AC", "#BA9FCC", "#9875B7", 
               "#764CA1", "#8B6899", "#C0AD99", "#F5F299", "#E8CE78", 
               "#CC9350", "#B15928")

# Create Sankey diagram
ggplot(dat, aes(axis1 = pog_I1.4, axis2 = pog_I2.0, axis3 = pog_I3.0, y = size)) +
  geom_alluvium(aes(fill = as.character(pog_I1.4)), width = 1/12, alpha = 0.7) +
  scale_fill_manual(values = my_colors, guide = "none") +
  geom_stratum(width = 1/12, color = "grey") +
  geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
  scale_x_continuous(
    breaks = 1:3, 
    labels = c(
      paste0("I=1.4\n", length(unique(dat$pog_I1.4)), ' POGs'), 
      paste0("I=2\n", length(unique(dat$pog_I2.0)), ' POGs'), 
      paste0("I=3.0\n", length(unique(dat$pog_I3.0)), ' POGs')
    )
  ) +
  labs(x = "Inflation Parameter", y = "Number of Proteins") +
  theme_minimal() +
  theme(
    legend.position = 'none', 
    panel.grid = element_blank(),
    
    # --- Key modifications ---
    # size = 12 can be adjusted as needed; larger numbers mean larger text
    # face = "bold" can be used for bold text
    axis.text.x = element_text(size = 12,color = "black",vjust = 7),
    axis.text.y = element_text(size = 12,color = "black"),
    # If you want to change the size of the x-axis title (Inflation Parameter) as well:
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14)
    # ------------------
  )

ggsave(file.path("outputs/Figure2", "Figure2C_Sankey_MCL_Inflation.png"), width=10, height=10, dpi=300)
ggsave(file.path("outputs/Figure2", "Figure2C_Sankey_MCL_Inflation.pdf"), width=10, height=10, dpi=300)
