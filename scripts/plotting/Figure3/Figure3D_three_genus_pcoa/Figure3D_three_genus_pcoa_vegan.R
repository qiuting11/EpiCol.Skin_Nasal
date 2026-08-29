#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(vegan)
  library(ggplot2)
})

`%||%` <- function(x, y) if (is.null(x)) y else x
script_path <- sys.frame(1)$ofile %||% commandArgs(trailingOnly = FALSE)[grep("^--file=", commandArgs(trailingOnly = FALSE))][1]
script_path <- sub("^--file=", "", script_path)
repo_root <- normalizePath(file.path(dirname(script_path), "..", "..", "..", ".."), mustWork = FALSE)
if (is.na(script_path) || !file.exists(file.path(repo_root, "README.md"))) {
  repo_root <- normalizePath(getwd(), mustWork = TRUE)
}

args <- commandArgs(trailingOnly = TRUE)
input_file <- if (length(args) >= 1) args[[1]] else file.path(repo_root, "data", "external", "figure3D_three_genus_pcoa", "mag_cf_long_with_core_phenotypes.tsv")
out_dir <- if (length(args) >= 2) args[[2]] else file.path(repo_root, "outputs", "Figure3", "Figure3D_three_genus_pcoa")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

n_perm <- as.integer(Sys.getenv("N_PERM", "999"))
if (is.na(n_perm) || n_perm < 1) stop("Invalid N_PERM")

target_genera <- c(
  "Staphylococcus" = "Staphylococcus",
  "Corynebacterium" = "Corynebacterium",
  "Cutibacterium_or_Propionibacterium" = "Cutibacterium"
)

plot_colors <- c(
  "Staphylococcus" = "#3B7EA1",
  "Corynebacterium" = "#D97742",
  "Cutibacterium" = "#4E9F6D"
)

log_msg <- function(...) {
  message(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", paste(..., collapse = ""))
  flush.console()
}

negative_eigen_summary <- function(eig) {
  pos <- eig[eig > 0]
  neg <- eig[eig < 0]
  total_abs <- sum(abs(eig))
  data.frame(
    n_eigenvalues = length(eig),
    n_positive = length(pos),
    n_negative = length(neg),
    positive_sum = sum(pos),
    negative_abs_sum = sum(abs(neg)),
    negative_abs_fraction = ifelse(total_abs == 0, NA_real_, sum(abs(neg)) / total_abs),
    min_eigenvalue = min(eig),
    stringsAsFactors = FALSE
  )
}

make_presence_matrix <- function(dat) {
  rows <- aggregate(presence_binary ~ MAG + CF_id, data = dat, FUN = max)
  xtabs(presence_binary ~ MAG + CF_id, data = rows)
}

make_copy_matrix <- function(dat) {
  rows <- aggregate(copy_number ~ MAG + CF_id, data = dat, FUN = sum)
  xtabs(copy_number ~ MAG + CF_id, data = rows)
}

run_one <- function(matrix, meta, metric, label, prefix) {
  log_msg("Building ", metric, " distance for ", label)
  if (metric == "jaccard") {
    dist_obj <- vegdist(matrix, method = "jaccard", binary = TRUE)
  } else if (metric == "bray") {
    dist_obj <- vegdist(log1p(matrix), method = "bray")
  } else {
    stop("Unsupported metric: ", metric)
  }

  log_msg("Running vegan::adonis2 for ", label, " with ", n_perm, " permutations")
  set.seed(20260827)
  ad <- adonis2(dist_obj ~ genus_display, data = meta, permutations = n_perm)

  log_msg("Running vegan::betadisper for ", label)
  bd <- betadisper(dist_obj, group = meta$genus_display, type = "centroid")

  log_msg("Running vegan::permutest.betadisper for ", label, " with ", n_perm, " permutations")
  set.seed(20260828)
  bd_perm <- permutest(bd, permutations = n_perm)

  eig <- bd$eig
  eig_pos <- eig[eig > 0]
  explained <- eig_pos / sum(eig_pos)
  coords <- as.data.frame(bd$vectors[, 1:2, drop = FALSE])
  names(coords) <- c("PCoA1", "PCoA2")

  coord_out <- cbind(meta, coords, distance_to_genus_centroid = bd$distances)
  write.table(coord_out, file.path(out_dir, paste0(prefix, "_coordinates.tsv")), sep = "\t", quote = FALSE, row.names = FALSE)

  eig_summary <- negative_eigen_summary(eig)
  eig_summary$analysis <- label
  eig_summary$metric <- metric
  write.table(eig_summary, file.path(out_dir, paste0(prefix, "_negative_eigen_summary.tsv")), sep = "\t", quote = FALSE, row.names = FALSE)

  stats <- data.frame(
    analysis = label,
    metric = metric,
    n_mags = nrow(meta),
    n_features = ncol(matrix),
    pcoa1_explained_positive_eig = explained[1],
    pcoa2_explained_positive_eig = explained[2],
    permanova_F = ad$F[1],
    permanova_R2 = ad$R2[1],
    permanova_p = ad$`Pr(>F)`[1],
    betadisper_F = bd_perm$tab$F[1],
    betadisper_p = bd_perm$tab$`Pr(>F)`[1],
    negative_eigenvalues = eig_summary$n_negative,
    negative_abs_fraction = eig_summary$negative_abs_fraction,
    permutations = n_perm,
    stringsAsFactors = FALSE
  )

  p <- ggplot(coord_out, aes(PCoA1, PCoA2, color = genus_display)) +
    geom_hline(yintercept = 0, color = "grey82", linewidth = 0.3) +
    geom_vline(xintercept = 0, color = "grey82", linewidth = 0.3) +
    geom_point(alpha = 0.55, size = 1.1) +
    scale_color_manual(values = plot_colors) +
    labs(
      title = label,
      x = sprintf("PCoA1 (%.1f%%)", explained[1] * 100),
      y = sprintf("PCoA2 (%.1f%%)", explained[2] * 100),
      color = NULL,
      caption = sprintf("PERMANOVA R2=%.3f, P<=%.3f; betadisper P<=%.3f; %d permutations", stats$permanova_R2, stats$permanova_p, stats$betadisper_p, n_perm)
    ) +
    theme_classic(base_size = 11) +
    theme(plot.title = element_text(face = "bold"), legend.position = c(0.82, 0.82), legend.background = element_rect(fill = "white", color = "grey85"))

  ggsave(file.path(out_dir, paste0(prefix, ".pdf")), p, width = 5.2, height = 4.6, units = "in")
  ggsave(file.path(out_dir, paste0(prefix, ".png")), p, width = 5.2, height = 4.6, units = "in", dpi = 300)
  stats
}

if (!file.exists(input_file)) {
  stop("Input file not found: ", input_file, "\nPlace the file under data/external/figure3D_three_genus_pcoa/ or pass input path as the first command-line argument.")
}

log_msg("Reading input: ", input_file)
dat <- read.delim(input_file, sep = "\t", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
required <- c("MAG", "CF_id", "copy_number", "target_genus_group", "Species")
missing <- setdiff(required, names(dat))
if (length(missing) > 0) stop("Missing required columns: ", paste(missing, collapse = ", "))

dat <- dat[dat$target_genus_group %in% names(target_genera), , drop = FALSE]
dat$genus_display <- unname(target_genera[dat$target_genus_group])
dat$copy_number <- as.numeric(dat$copy_number)
dat$presence_binary <- as.integer(dat$copy_number > 0)

presence <- make_presence_matrix(dat)
copy_number <- make_copy_matrix(dat)
common_mags <- intersect(rownames(presence), rownames(copy_number))
presence <- presence[common_mags, , drop = FALSE]
copy_number <- copy_number[common_mags, , drop = FALSE]

meta_cols <- c("MAG", "genus_display", "target_genus_group", "Species")
if ("source_group" %in% names(dat)) meta_cols <- c(meta_cols, "source_group")
meta <- dat[meta_cols]
meta <- meta[!duplicated(meta$MAG), , drop = FALSE]
meta <- meta[match(common_mags, meta$MAG), , drop = FALSE]
stopifnot(identical(meta$MAG, rownames(presence)))

sample_counts <- as.data.frame(table(meta$genus_display), stringsAsFactors = FALSE)
names(sample_counts) <- c("genus_display", "n_mags")
write.table(sample_counts, file.path(out_dir, "Figure3D_three_genus_pcoa_sample_counts.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

stats <- rbind(
  run_one(presence, meta, "jaccard", "A. eCF presence/absence profiles", "Figure3D_three_genus_presence_jaccard_pcoa"),
  run_one(copy_number, meta, "bray", "B. eCF copy-number profiles", "Figure3D_three_genus_copy_braycurtis_pcoa")
)

write.table(stats, file.path(out_dir, "Figure3D_three_genus_pcoa_statistics.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
print(stats)
log_msg("Wrote outputs to: ", normalizePath(out_dir, winslash = "/"))
