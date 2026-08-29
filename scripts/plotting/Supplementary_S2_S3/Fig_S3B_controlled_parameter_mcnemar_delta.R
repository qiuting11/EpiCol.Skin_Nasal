#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(cowplot)
  library(ggplot2)
})

repo_root <- Sys.getenv("PAPER_SUBMISSION_CODE", unset = normalizePath(getwd(), winslash = "/", mustWork = TRUE))
root <- file.path(repo_root, "data", "external", "supplementary_s2_s3", "human_hq_mag_cf_matrix")
out_dir <- file.path(repo_root, "outputs", "Supplementary_S2_S3")
plot_table_dir <- file.path(out_dir, "plot_tables")
tables_dir <- file.path(root, "controlled_parameter_comparison", "tables")
presence_tables_dir <- file.path(root, "main_strict_sensitive_result", "cf_presence_absence_analysis", "tables")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_table_dir, recursive = TRUE, showWarnings = FALSE)

source(file.path(repo_root, "scripts", "plotting", "Figure3", "Figure3_color_config.R"))

read_tsv_local <- function(path) {
  read.delim(path, sep = "\t", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
}

save_plot <- function(plot, path_base, width, height) {
  ggsave(paste0(path_base, ".pdf"), plot, width = width, height = height, units = "in", device = cairo_pdf)
  ggsave(paste0(path_base, ".png"), plot, width = width, height = height, units = "in", dpi = 300)
}

method_files <- list(
  strict_fast = file.path(root, "sensitivity_strict_fast", "outputs", "mag_cf_presence_with_phenotypes.tsv"),
  strict_sensitive = file.path(root, "sensitivity_strict_sensitive", "outputs", "mag_cf_presence_with_phenotypes.tsv"),
  broad = file.path(root, "sensitivity_paper_like", "outputs", "mag_cf_presence_with_phenotypes.tsv")
)

presence <- lapply(method_files, read_tsv_local)
cf_cols <- grep("^CF", names(presence$strict_fast), value = TRUE)

cf_reference <- read_tsv_local(file.path(presence_tables_dir, "cf_presence_taxonomy_summary_for_plot.tsv"))
cf_reference <- cf_reference[!duplicated(cf_reference$CF_id), c("CF_id", "eCM", "eCM_num")]
cf_reference$cf_major <- as.integer(sub("^CF([0-9]+).*$", "\\1", cf_reference$CF_id))
cf_reference$cf_minor <- -1L
has_minor <- grepl("_", cf_reference$CF_id)
cf_reference$cf_minor[has_minor] <- as.integer(sub("^CF[0-9]+_([0-9]+)$", "\\1", cf_reference$CF_id[has_minor]))
cf_reference <- cf_reference[order(cf_reference$eCM_num, cf_reference$cf_major, cf_reference$cf_minor, cf_reference$CF_id), ]
cf_order <- cf_reference$CF_id
extra_cf <- setdiff(cf_cols, cf_order)
if (length(extra_cf) > 0) {
  cf_order <- c(cf_order, sort(extra_cf))
}
cf_labels <- setNames(sub("^CF", "eCF", cf_order), cf_order)

module_bar <- cf_reference
module_bar$CF_id <- factor(module_bar$CF_id, levels = cf_order)
module_bar$eCM_label <- sub("^module", "eCM", module_bar$eCM)
module_bar$module_fill <- unname(full_module_colors[module_bar$eCM])
module_bar$module_fill[is.na(module_bar$module_fill)] <- "#D9D9D9"

genus_label_map <- c(
  "Staphylococcus" = "Staphylococcus",
  "Corynebacterium" = "Corynebacterium",
  "Cutibacterium_or_Propionibacterium" = "Cutibacterium",
  "Cutibacterium" = "Cutibacterium"
)
genus_order <- c("Cutibacterium", "Corynebacterium", "Staphylococcus")

make_long <- function(df, method_id) {
  keep <- c("MAG", "target_genus_group", cf_cols)
  df <- df[, keep]
  out <- reshape(
    df,
    varying = cf_cols,
    v.names = "presence",
    timevar = "CF_id",
    times = cf_cols,
    idvar = c("MAG", "target_genus_group"),
    direction = "long"
  )
  rownames(out) <- NULL
  out$method_id <- method_id
  out$presence <- as.integer(out$presence > 0)
  out
}

long <- rbind(
  make_long(presence$strict_fast, "strict_fast"),
  make_long(presence$strict_sensitive, "strict_sensitive"),
  make_long(presence$broad, "broad")
)

calc_contrast <- function(method_a, method_b, contrast_label) {
  a <- long[long$method_id == method_a, c("MAG", "target_genus_group", "CF_id", "presence")]
  b <- long[long$method_id == method_b, c("MAG", "CF_id", "presence")]
  names(a)[names(a) == "presence"] <- "presence_a"
  names(b)[names(b) == "presence"] <- "presence_b"
  merged <- merge(a, b, by = c("MAG", "CF_id"), all = FALSE)
  pieces <- by(merged, list(merged$target_genus_group, merged$CF_id), function(x) {
    a_only <- sum(x$presence_a == 1 & x$presence_b == 0, na.rm = TRUE)
    b_only <- sum(x$presence_a == 0 & x$presence_b == 1, na.rm = TRUE)
    discordant <- a_only + b_only
    p_value <- if (discordant == 0) 1 else binom.test(min(a_only, b_only), discordant, p = 0.5, alternative = "two.sided")$p.value
    data.frame(
      contrast = contrast_label,
      method_a = method_a,
      method_b = method_b,
      target_genus_group = x$target_genus_group[1],
      CF_id = x$CF_id[1],
      n_mags = length(unique(x$MAG)),
      prevalence_a = mean(x$presence_a, na.rm = TRUE),
      prevalence_b = mean(x$presence_b, na.rm = TRUE),
      delta_prevalence = mean(x$presence_b, na.rm = TRUE) - mean(x$presence_a, na.rm = TRUE),
      a_present_b_absent = a_only,
      a_absent_b_present = b_only,
      discordant_pairs = discordant,
      p_value = p_value
    )
  })
  do.call(rbind, pieces)
}

mcnemar <- rbind(
  calc_contrast("strict_fast", "strict_sensitive", "Strict-sensitive - Strict-fast"),
  calc_contrast("strict_sensitive", "broad", "Broad homolog-recovery - Strict-sensitive")
)

mcnemar$FDR <- ave(mcnemar$p_value, mcnemar$contrast, FUN = function(x) p.adjust(x, method = "BH"))
mcnemar$significant <- !is.na(mcnemar$FDR) & mcnemar$FDR < 0.05
mcnemar$significance_mark <- ifelse(mcnemar$significant, "*", "")
mcnemar$target_genus_display <- genus_label_map[mcnemar$target_genus_group]
mcnemar$target_genus_display[is.na(mcnemar$target_genus_display)] <- mcnemar$target_genus_group[is.na(mcnemar$target_genus_display)]
mcnemar$target_genus_display <- factor(mcnemar$target_genus_display, levels = genus_order)
mcnemar$CF_id <- factor(mcnemar$CF_id, levels = cf_order)
mcnemar$contrast <- factor(mcnemar$contrast, levels = c("Strict-sensitive - Strict-fast", "Broad homolog-recovery - Strict-sensitive"))

write.table(mcnemar, file.path(plot_table_dir, "controlled_genus_paired_mcnemar_presence_delta.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

make_module_strip <- function() {
  ggplot(module_bar, aes(x = CF_id, y = 1, fill = module_fill)) +
    geom_tile(color = "white", linewidth = 0.25) +
    geom_text(aes(label = eCM_label), size = 2, color = "black") +
    scale_fill_identity() +
    scale_x_discrete(limits = cf_order, labels = cf_labels, drop = FALSE, expand = expansion(add = 0.6)) +
    scale_y_continuous(expand = c(0, 0)) +
    labs(x = NULL, y = NULL) +
    theme_void(base_size = 8) +
    theme(plot.margin = margin(0, 50, 0, 72))
}

theme_delta <- theme_bw(base_size = 9) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7, color = "black"),
    axis.text.y = element_text(size = 8, color = "black"),
    strip.background = element_rect(fill = "grey92", color = "grey70"),
    strip.text.y = element_text(face = "bold", size = 8),
    legend.position = "bottom",
    legend.box = "horizontal",
    plot.title = element_blank(),
    plot.subtitle = element_blank()
  )

p_delta <- ggplot(mcnemar, aes(x = CF_id, y = target_genus_display, fill = delta_prevalence)) +
  geom_tile(color = "white", linewidth = 0.25) +
  geom_text(aes(label = significance_mark), size = 3.2, color = "black") +
  facet_grid(contrast ~ ., switch = "y") +
  scale_fill_gradient2(
    low = "#2166AC",
    mid = "white",
    high = "#B2182B",
    midpoint = 0,
    limits = c(-1, 1),
    name = "Delta prevalence\nmethod B - method A",
    labels = scales::percent_format(accuracy = 1)
  ) +
  scale_x_discrete(limits = cf_order, labels = cf_labels, drop = FALSE, expand = expansion(add = 0.6)) +
  labs(x = "eCF", y = "Genus") +
  theme_delta

title_plot <- cowplot::ggdraw() +
  cowplot::draw_label("Paired changes in genus-stratified eCF presence under controlled DIAMOND parameter sets", x = 0, y = 0.68, hjust = 0, fontface = "bold", size = 11) +
  cowplot::draw_label("Fill shows paired prevalence change; * marks McNemar FDR < 0.05 within each contrast.", x = 0, y = 0.25, hjust = 0, size = 8, color = "grey30") +
  theme(plot.margin = margin(2, 48, 0, 72))

p_out <- cowplot::plot_grid(title_plot, make_module_strip(), p_delta, ncol = 1, rel_heights = c(0.18, 0.10, 1), align = "v", axis = "lr")
save_plot(p_out, file.path(out_dir, "Fig_S3B_controlled_genus_paired_mcnemar_delta_prevalence_heatmap"), 12.2, 5.8)

cat("McNemar table written:", file.path(plot_table_dir, "controlled_genus_paired_mcnemar_presence_delta.tsv"), "\n")
cat("McNemar delta figure written:", file.path(out_dir, "Fig_S3B_controlled_genus_paired_mcnemar_delta_prevalence_heatmap"), "\n")
cat("Significant cells by contrast:\n")
print(aggregate(significant ~ contrast, data = mcnemar, FUN = sum))
