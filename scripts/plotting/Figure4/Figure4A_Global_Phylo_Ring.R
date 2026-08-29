library(ape)
library(ggplot2)
library(ggtree)
library(ggtreeExtra)
library(ggnewscale)
library(dplyr)
library(readr)
library(stringr)
library(tibble)

# -----------------------------
# Config
# -----------------------------
input_tree <- "data/figure4/gtdb_r226_phylum_reduced.tree"
input_phylum_size <- "data/figure4/phylum_size.csv"
input_cf_coverage <- "data/figure4/cf_phylum_coverage_long.csv"
input_modules <- "data/figure4/top0.2_I3.0_modules.txt"
output_dir <- "outputs/Figure4"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

out_pdf <- file.path(output_dir, "Figure4A_Global_Phylo_Ring.pdf")
out_png <- file.path(output_dir, "Figure4A_Global_Phylo_Ring.png")
out_summary <- file.path(output_dir, "Figure4A_phylum_summary.tsv")
out_missing_cf <- file.path(output_dir, "Figure4A_missing_cf_list.tsv")

# -----------------------------
# Helpers
# -----------------------------
norm_phylum <- function(x) {
  x <- trimws(as.character(x))
  ifelse(str_detect(x, "^p__"), x, paste0("p__", x))
}

# -----------------------------
# Load inputs
# -----------------------------
stopifnot(file.exists(input_tree), file.exists(input_phylum_size), file.exists(input_cf_coverage), file.exists(input_modules))

tree <- read.tree(input_tree)

phylum_size_raw <- read.csv(input_phylum_size, check.names = FALSE, stringsAsFactors = FALSE)
if ("phylum" %in% names(phylum_size_raw)) {
  phylum_size <- phylum_size_raw %>%
    transmute(phylum = norm_phylum(phylum), size = as.numeric(size))
} else {
  phylum_size <- tibble(
    phylum = norm_phylum(phylum_size_raw[[1]]),
    size = as.numeric(phylum_size_raw[[2]])
  )
}

cf_coverage <- read.csv(input_cf_coverage, stringsAsFactors = FALSE) %>%
  mutate(
    phylum = norm_phylum(phylum),
    CF_id = as.character(CF_id),
    percentage = as.numeric(percentage),
    positive_species = as.numeric(positive_species),
    total_species_in_phylum = as.numeric(total_species_in_phylum)
  )

modules <- read_delim(input_modules, delim = "\t", show_col_types = FALSE) %>%
  select(module, cf_id) %>%
  distinct() %>%
  mutate(mod_num = as.numeric(str_extract(module, "\\d+")))

# Core fix: remove CFs that did not appear in the search results
missing_cfs_known <- c("CF15_19", "CF17_23", "CF9_12")
cf_order <- modules %>%
  filter(!cf_id %in% missing_cfs_known) %>%
  arrange(mod_num, cf_id) %>%
  pull(cf_id) %>%
  unique()

cf_coverage <- cf_coverage %>%
  filter(CF_id %in% cf_order) %>%
  mutate(CF_id = factor(CF_id, levels = cf_order))

# -----------------------------
# Summary stats
# -----------------------------
mean_cf_count_refined <- cf_coverage %>%
  group_by(phylum) %>%
  summarize(
    total_hits = sum(positive_species, na.rm = TRUE),
    total_sp = first(total_species_in_phylum),
    min_capable_sp = max(positive_species, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    mean_count_active = if_else(min_capable_sp > 0, total_hits / min_capable_sp, 0),
    rank_val = if_else(total_sp > 5, rank(-mean_count_active, ties.method = "first"), NA_real_),
    is_top20 = if_else(!is.na(rank_val) & rank_val <= 20, "Highlight", "Normal")
  )

phylum_size_data <- phylum_size

summary_tbl <- phylum_size %>%
  left_join(mean_cf_count_refined, by = "phylum") %>%
  arrange(desc(size))
write.table(summary_tbl, out_summary, sep = "\t", row.names = FALSE, quote = FALSE)

missing_cf <- setdiff(sort(unique(modules$cf_id)), sort(unique(as.character(cf_coverage$CF_id))))
write.table(
  data.frame(missing_cf = missing_cf),
  out_missing_cf,
  sep = "\t", row.names = FALSE, quote = FALSE
)

# -----------------------------
# Plot Figure4A
# -----------------------------
p <- ggtree(tree, layout = "fan", open.angle = 5, linewidth = 0.2) +
  geom_hline(yintercept = seq_len(length(tree$tip.label)), linetype = "dashed", color = "grey90", linewidth = 0.1)

p <- p %<+% phylum_size_data +
  geom_tippoint(aes(size = size), color = "grey30", alpha = 0.45) +
  geom_tiplab(
    aes(label = paste0(label, " (", size, ")"), color = (size > 10000)),
    align = TRUE, linetype = "dotted", size = 1.5, offset = 5
  ) +
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "black"), guide = "none") +
  scale_size_continuous(range = c(0.5, 4), name = "Species Total")

p <- p +
  geom_fruit(
    data = cf_coverage,
    geom = geom_tile,
    mapping = aes(y = phylum, x = CF_id, fill = percentage),
    offset = 0.08, pwidth = 1.2, color = "white", linewidth = 0.05,
    axis.params = list(axis = "x", text.angle = 90, text.size = 1.2, vjust = 0.5, text = function(x) sub("^CF", "eCF", x)),
    grid.params = list(linetype = "dotted", color = "grey80", linewidth = 0.1)
  ) +
  scale_fill_viridis_c(option = "D", name = "Coverage (%)") +
  new_scale_fill()

p <- p +
  geom_fruit(
    data = mean_cf_count_refined,
    geom = geom_bar,
    mapping = aes(y = phylum, x = mean_count_active, fill = is_top20),
    stat = "identity", orientation = "y",
    offset = 0.15, pwidth = 0.4,
    axis.params = list(axis = "x", title = "Mean eCFs", title.size = 2, text.size = 1.5)
  ) +
  scale_fill_manual(values = c("Highlight" = "#E41A1C", "Normal" = "grey70"), name = "Status") +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 7)
  )

suppressWarnings(ggsave(out_pdf, plot = p, width = 16, height = 16, device = cairo_pdf))
suppressWarnings(ggsave(out_png, plot = p, width = 20, height = 20, dpi = 300))

message("Figure4A done: ", out_pdf)
message("Figure4A done: ", out_png)
message("Summary: ", out_summary)
message("Missing eCF list: ", out_missing_cf)
