library(dplyr)
library(readr)
library(tidyr)
library(ggplot2)
library(stringr)
library(forcats)

# -----------------------------
# Config
# -----------------------------
input_cf_coverage_pos <- "E:/data/mmseqs/sp/04_Final_Metrics/cf_phylum_coverage_long.csv"   # positive-only long table
input_phylum_size <- "E:/data/mmseqs/sp/04_Final_Metrics/phylum_size.csv"
input_modules <- "E:/data/gephe_output_11/analysis/output_test_1/top0.2_I3.0_modules.txt"
output_dir <- "E:/data/R/gephe_R/results/figures/Figure4_Final"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

out_full_long <- file.path(output_dir, "Figure4B_plotdata_full_withZero.tsv")
out_pos_long  <- file.path(output_dir, "Figure4B_plotdata_positive.tsv")
out_matrix    <- file.path(output_dir, "Figure4B_phylum_cf_matrix_withZero.tsv")
out_summary   <- file.path(output_dir, "Figure4B_coverage_summary.tsv")

out_full_pdf <- file.path(output_dir, "Figure4B_fullCoverage_withZero.pdf")
out_full_png <- file.path(output_dir, "Figure4B_fullCoverage_withZero.png")
out_pos_pdf  <- file.path(output_dir, "Figure4B_positiveOnly.pdf")
out_pos_png  <- file.path(output_dir, "Figure4B_positiveOnly.png")

norm_phylum <- function(x){
  x <- trimws(as.character(x))
  ifelse(str_detect(x, "^p__"), x, paste0("p__", x))
}

# -----------------------------
# Load
# -----------------------------
pos <- read.csv(input_cf_coverage_pos, stringsAsFactors = FALSE) %>%
  mutate(
    phylum = norm_phylum(phylum),
    CF_id = as.character(CF_id),
    percentage = as.numeric(percentage)
  )

mods <- read_delim(input_modules, delim = "\t", show_col_types = FALSE) %>%
  select(module, cf_id) %>% distinct() %>%
  mutate(mod_num = as.numeric(str_extract(module, "\\d+")))

# Core fix: remove CFs that did not appear in the search results
missing_cfs_known <- c("CF15_19", "CF17_23", "CF9_12")
cf_order <- mods %>% 
  filter(!cf_id %in% missing_cfs_known) %>%
  arrange(mod_num, cf_id) %>% 
  pull(cf_id) %>% 
  unique()

ph_raw <- read.csv(input_phylum_size, check.names = FALSE, stringsAsFactors = FALSE)
if ("phylum" %in% names(ph_raw)) {
  ph <- ph_raw %>% transmute(phylum = norm_phylum(phylum), total_sp = as.numeric(size))
} else {
  ph <- tibble(phylum = norm_phylum(ph_raw[[1]]), total_sp = as.numeric(ph_raw[[2]]))
}

# -----------------------------
# Rebuild full matrix with zeros
# -----------------------------
all_grid <- tidyr::expand_grid(phylum = unique(ph$phylum), CF_id = unique(cf_order))

full <- all_grid %>%
  left_join(ph, by = "phylum") %>%
  left_join(pos %>% select(phylum, CF_id, percentage, positive_species, total_species_in_phylum),
            by = c("phylum","CF_id")) %>%
  mutate(
    percentage = if_else(is.na(percentage), 0, percentage),
    total_species_in_phylum = if_else(is.na(total_species_in_phylum), total_sp, total_species_in_phylum),
    positive_species = if_else(is.na(positive_species), round(percentage * total_species_in_phylum / 100), positive_species)
  ) %>%
  left_join(mods, by = c("CF_id" = "cf_id")) %>%
  mutate(
    mod_num = as.numeric(str_extract(module, "\\d+")),
    CF_id = factor(CF_id, levels = cf_order),
    module = fct_reorder(module, mod_num, .na_rm = FALSE)
  )

pos_only <- full %>% filter(percentage > 0)

# wide matrix output
mat <- full %>% select(phylum, CF_id, percentage) %>%
  mutate(CF_id = as.character(CF_id)) %>%
  tidyr::pivot_wider(names_from = CF_id, values_from = percentage, values_fill = 0)

# summary
sum_tbl <- full %>%
  group_by(CF_id, module) %>%
  summarise(
    n_phylum = n(),
    n_zero = sum(percentage == 0),
    frac_zero = mean(percentage == 0),
    median_cov = median(percentage),
    iqr_cov = IQR(percentage),
    .groups = "drop"
  ) %>%
  arrange(desc(median_cov))

write.table(full, out_full_long, sep = "\t", row.names = FALSE, quote = FALSE)
write.table(pos_only, out_pos_long, sep = "\t", row.names = FALSE, quote = FALSE)
write.table(mat, out_matrix, sep = "\t", row.names = FALSE, quote = FALSE)
write.table(sum_tbl, out_summary, sep = "\t", row.names = FALSE, quote = FALSE)

# -----------------------------
# Plot
# -----------------------------
my_colors <- c(
  "#72190E", "#f94144", "#F7557F", "#7C417F", "#A05992", "#f3722c",
  "#f8961e", "#f9c74f", "#F7F056", "#C2B923", "#A9D88C", "#90be6d",
  "#43aa8b", "#277da1", "#1BA3C6", "#B9F2F0", "#4F7CBA", "#5289C7"
)
mods_u <- sort(unique(as.character(full$module)))
if(length(mods_u) > length(my_colors)) my_colors <- rep(my_colors, length.out = length(mods_u))
names(my_colors) <- mods_u

mk <- function(df, title, subtitle){
  ggplot(df, aes(x = CF_id, y = percentage)) +
    geom_violin(aes(fill = module), alpha = 0.2, color = NA, scale = "width") +
    geom_boxplot(width = 0.18, outlier.shape = NA, alpha = 0.35, color = "grey30") +
    geom_jitter(aes(size = total_sp, color = module), width = 0.18, alpha = 0.30) +
    scale_fill_manual(values = my_colors, name = "CM(module)") +
    scale_color_manual(values = my_colors, name = "CM(module)") +
    scale_size_continuous(range = c(0.8, 4.8), name = "Phylum size") +
    labs(title = title, subtitle = subtitle, x = "CF", y = "Species coverage in phylum (%)") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8),
          panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
          legend.position = "right")
}

p_full <- mk(full, "Figure4B (MMseqs): Full coverage distribution", "Includes zero-coverage phylum?CF combinations")
p_pos <- mk(pos_only, "Figure4B (MMseqs): Positive-only distribution", "Only combinations with coverage > 0")

ggsave(out_full_pdf, p_full, width = 15, height = 7, device = cairo_pdf)
ggsave(out_full_png, p_full, width = 15, height = 7, dpi = 300)
ggsave(out_pos_pdf, p_pos, width = 15, height = 7, device = cairo_pdf)
ggsave(out_pos_png, p_pos, width = 15, height = 7, dpi = 300)

message("Done Figure4B rebuild with zeros")
message(out_full_long)
message(out_summary)
