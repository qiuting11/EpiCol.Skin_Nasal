library(tidyverse)  

# Read CSV file  
result_df <- read_csv('E:/data/gephe_output_11/association/result_11_1230/mi_zscore.csv')  
breaks <- c(0.99, 0.995, 0.996, 0.997, 0.998, 0.999)  
quantile_values <- quantile(result_df$mi_z, breaks, na.rm = TRUE)  

# Create label data frame  
breaks_df <- data.frame(  
  breaks = quantile_values,  
  y = Inf,  
  label = paste0(round((1-breaks)*100, 2), '%')  
)  

# Expand x-axis range to ensure all thresholds are visible  
x_min <- min(result_df$mi_z, na.rm = TRUE)  
# 1. Determine global maximum to avoid hard truncation
x_max_real <- max(result_df$mi_z, na.rm = TRUE)

# 2. Improved plotting function
plot_mi_hist_fixed <- function(dat, breaks_data, title = '') {
  # Get the maximum frequency of the current data for dynamic label height setting
  # This avoids warnings caused by using Inf
  max_count <- max(table(cut(dat$mi_z, breaks = 100)))
  
  ggplot(dat, aes(mi_z)) + 
    geom_histogram(bins = 100, fill =  "#000000", alpha = 0.8) +
    # Use coord_cartesian instead of xlim, which preserves data and only performs window scaling
    coord_cartesian(xlim = c(x_min, x_max_real)) +
    geom_vline(xintercept = breaks_data$breaks, 
               color = "red", linetype = "dashed", linewidth = 0.6) +
    # Use check_overlap from the text package or adjust position
    geom_text(aes(x = breaks, y = max_count, label = label), 
              data = breaks_data,
              vjust = -0.5, # Offset upwards to avoid overlapping the scale lines
              angle = 90,   # Rotate 90 degrees, suitable for dense threshold labeling
              size = 3.5, 
              color = "red") +
    labs(x = "Mutual Information Z-score", 
         y = 'Count',
         title = title) +
    theme_bw() +
    theme(panel.grid.minor = element_blank()) # Cleaner theme
}

# 3. Execute plotting
mi_z_fixed <- plot_mi_hist_fixed(result_df, breaks_df, title = "MI Z-score Distribution")
print(mi_z_fixed)
ggsave("FigureS2B_MI_Zscore_Distribution_FullQuantiles.png", mi_z_fixed, dpi = 300, width = 10, height = 6)

# 1. Re-filter to show only the 0.998 quantile
target_break <- 0.998
q_val_0998 <- quantile(result_df$mi_z, target_break, na.rm = TRUE)

# 2. Create a label data frame containing only a single line
# Here label contains the percentage and the specific Z-score value
single_break_df <- data.frame(
  breaks = q_val_0998,
  y = Inf, 
  label = paste0("Top 0.2%\n (Z >= ", round(q_val_0998, 3), ")")
)

# 3. Update plotting function (optimized for single line)
plot_mi_single_line <- function(dat, break_data, title = '') {
  # Automatically get the maximum y-axis value to determine text position
  # Calculate silently to avoid complex calculations on tens of millions of data points
  max_y <- max(ggplot_build(
    ggplot(dat, aes(mi_z)) + geom_histogram(bins = 100)
  )$data[[1]]$count)
  
  ggplot(dat, aes(mi_z)) + 
    geom_histogram(bins = 100, fill = "#000000", alpha = 0.8) +
    coord_cartesian(xlim = c(x_min, x_max_real)) +
    # Draw the 0.998 line
    geom_vline(xintercept = break_data$breaks, 
               color = "red", linetype = "dashed", linewidth = 0.8) +
    # Text annotation: displayed near the highest point
    geom_text(aes(x = breaks, y = max_y * 0.9, label = label), 
              data = break_data,
              hjust = -0.1,  # Offset slightly to the right of the line
              size = 5, 
              fontface = "bold",
              color = "red") +
    labs(x = "Mutual Information Z-score", 
         y = 'Count',
         title = title) +
    theme_bw()
}

# Execute plotting
mi_z_0998_only <- plot_mi_single_line(result_df, single_break_df, title = "MI Z-score Distribution (Highlight 0.998)")
print(mi_z_0998_only)
ggsave("FigureS2B_MI_Zscore_Distribution_Alt.png", mi_z_0998_only, dpi = 300, width = 8, height = 6)
ggsave("FigureS2B_MI_Zscore_Distribution_Alt.pdf", mi_z_0998_only, dpi = 300, width = 8, height = 6)

# 1. Calculate all high quantiles (keep the calculation process for future comparison)
top_breaks <- c(0.995, 0.996, 0.997, 0.998, 0.999)
top_quantiles <- quantile(result_df$mi_z, top_breaks, na.rm = TRUE)

# 2. Extract only 0.998 data and construct detailed labels
# This will display text like "0.2% (Z=12.45)" on the chart
target_idx <- which(top_breaks == 0.998)
single_break_df <- data.frame(
  breaks = top_quantiles[target_idx],
  y = Inf,
  label = paste0("Top 0.2%\n(Z=", round(top_quantiles[target_idx], 2), ")")
)

# 3. Filter Top 1% data for plotting (close-up view)
top_threshold_99 <- quantile(result_df$mi_z, 0.99, na.rm = TRUE)
top_data <- result_df[result_df$mi_z > top_threshold_99, ]

# 4. Plotting: show only the 0.998 line
g_top <- ggplot(top_data, aes(mi_z)) +   
  geom_histogram(bins = 100, fill ="#000000", alpha = 0.6) +  # Increase bins for finer detail in close-up
  # Draw the 0.998 vertical line
  geom_vline(xintercept = single_break_df$breaks, 
             color =  "red", linetype = "dashed", linewidth = 1.2) +  
  # Label percentage and Z-score value
  geom_text(aes(x = breaks, y = y, label = label), 
            data = single_break_df,   
            vjust = 1.5,      # Adjust vertical position to place it below the top
            hjust = -0.1,     # Offset to the right of the line to avoid overlapping
            size = 5, 
            fontface = "bold",
            color =  "red") +  
  labs(x = "Mutual Information Z-score", 
       y = 'Count',
       title = "MI Z-score Distribution (Top 1% Close-up)") +
  theme_bw() +
  # Automatically adjust X-axis range to ensure line and text are included
  coord_cartesian(clip = "off") 

# Save plots
ggsave("Figure2B_MI_Zscore_Distribution.png", g_top, dpi = 300, width = 8, height = 6)
ggsave("Figure2B_MI_Zscore_Distribution.pdf", g_top, dpi = 300, width = 8, height = 6)
