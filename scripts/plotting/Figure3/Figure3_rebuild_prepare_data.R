suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
})

script_dir <- "scripts/plotting/Figure3"
source(file.path(script_dir, "Figure3_color_config.R"))

paths <- list(
  modules = "data/figure3/top0.2_I3.0_modules.txt",
  eggnog = "data/external/eggnog_output_with_pog_multi_level.tsv",
  metadata = "data/figure3/metadata_for_analysis.tsv",
  pog_copy = "data/figure3/pog_copy_number_profiles.tsv",
  out_dir = "outputs/Figure3"
)
dir.create(paths$out_dir, recursive = TRUE, showWarnings = FALSE)

read_tsv <- function(path) {
  read.delim(path, sep = "\t", header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
}

normalize_habitat <- function(source, habitat = NA_character_) {
  habitat <- as.character(habitat)
  source <- as.character(source)
  dplyr::case_when(
    habitat %in% c("Skin", "Nasal", "Environment") ~ habitat,
    habitat %in% c("Nose", "nose", "nasal") ~ "Nasal",
    habitat %in% c("Environmental") ~ "Environment",
    source %in% c("HSMG", "UHSG") ~ "Skin",
    source %in% c("nose", "Nose", "nasal", "Nasal") ~ "Nasal",
    source %in% c("Aquatic", "Terrestrial", "other") ~ "Environment",
    TRUE ~ "Unknown"
  )
}

collapse_unique <- function(x, max_items = 4) {
  x <- unique(x[!is.na(x) & x != "" & x != "-"])
  if (length(x) == 0) return("")
  paste(head(x, max_items), collapse = "; ")
}

modules_df <- read_tsv(paths$modules) %>%
  mutate(
    module = factor(module, levels = module_order),
    functional_theme = unname(functional_theme_map[as.character(module)]),
    functional_theme_short = unname(functional_theme_short[functional_theme]),
    ecological_relevance = unname(functional_theme_relevance[functional_theme]),
    pog_row_name = paste0("POG_", as.numeric(gsub("pog", "", original_pog)))
  ) %>%
  arrange(module, cf_id)

all_cf_ids <- modules_df$cf_id

eggnog_df <- read_tsv(paths$eggnog) %>%
  filter(CF_id %in% all_cf_ids)

cf_annotation <- eggnog_df %>%
  group_by(CF_id) %>%
  summarise(
    preferred_name = collapse_unique(Preferred_name),
    description = collapse_unique(Description),
    kegg_ko = collapse_unique(KEGG_ko),
    cog_category = collapse_unique(COG_category),
    .groups = "drop"
  )

cf_info <- modules_df %>%
  left_join(cf_annotation, by = c("cf_id" = "CF_id")) %>%
  mutate(
    preferred_name = ifelse(is.na(preferred_name) | preferred_name == "", keyword, preferred_name),
    description = ifelse(is.na(description), "", description),
    kegg_ko = ifelse(is.na(kegg_ko), "", kegg_ko),
    cog_category = ifelse(is.na(cog_category), "", cog_category)
  )

figure3a_cf_table <- cf_info %>%
  transmute(
    functional_theme,
    functional_theme_short,
    ecological_relevance,
    module = as.character(module),
    CF = cf_id,
    original_pog,
    keyword,
    preferred_name,
    description,
    kegg_ko,
    cog_category,
    N,
    mi,
    mi_z,
    spearman_r
  )

figure3a_module_summary <- figure3a_cf_table %>%
  group_by(functional_theme, functional_theme_short, ecological_relevance, module) %>%
  summarise(
    n_cf = n_distinct(CF),
    representative_eCFs = paste(CF, collapse = "; "),
    representative_keywords = collapse_unique(keyword, max_items = 3),
    .groups = "drop"
  ) %>%
  arrange(match(module, module_order))

metadata_df <- read_tsv(paths$metadata) %>%
  mutate(
    Habitat = normalize_habitat(source, if ("Habitat" %in% names(.)) Habitat else NA_character_),
    Genus = if ("Genus" %in% names(.)) Genus else str_extract(classification, "(?<=g__)[^;]+"),
    Species = if ("Species" %in% names(.)) Species else str_extract(classification, "(?<=s__)[^;]*"),
    Genus = ifelse(is.na(Genus) | Genus == "", str_extract(classification, "(?<=g__)[^;]+"), Genus),
    Species = ifelse(is.na(Species) | Species == "", paste0(Genus, " sp."), Species)
  ) %>%
  filter(!is.na(Genus), Genus != "", !grepl("unclassified|unknown|Incertae_Sedis", Genus, ignore.case = TRUE))

pog_copy_mat <- as.matrix(read.table(paths$pog_copy, sep = "\t", header = TRUE, row.names = 1, check.names = FALSE))
storage.mode(pog_copy_mat) <- "numeric"

cf_copy_sum_mat <- matrix(0, nrow = length(all_cf_ids), ncol = ncol(pog_copy_mat))
rownames(cf_copy_sum_mat) <- all_cf_ids
colnames(cf_copy_sum_mat) <- colnames(pog_copy_mat)

for (cf in all_cf_ids) {
  pogs <- cf_info$pog_row_name[cf_info$cf_id == cf]
  pogs <- pogs[pogs %in% rownames(pog_copy_mat)]
  if (length(pogs) > 0) {
    cf_copy_sum_mat[cf, ] <- colSums(pog_copy_mat[pogs, , drop = FALSE], na.rm = TRUE)
  }
}

cf_copy_wide <- as.data.frame(t(cf_copy_sum_mat), check.names = FALSE) %>%
  tibble::rownames_to_column("taxon_oid") %>%
  inner_join(metadata_df, by = "taxon_oid")

human_surface_wide <- cf_copy_wide %>%
  filter(Habitat %in% c("Skin", "Nasal"))

genus_habitat_counts <- human_surface_wide %>%
  group_by(Genus) %>%
  summarise(
    n_skin = sum(Habitat == "Skin"),
    n_nasal = sum(Habitat == "Nasal"),
    n_genomes = n(),
    n_species = n_distinct(Species),
    .groups = "drop"
  ) %>%
  mutate(
    HabitatGroup = case_when(
      n_skin > 0 & n_nasal > 0 ~ "Shared",
      n_skin > 0 & n_nasal == 0 ~ "Skin-associated",
      n_nasal > 0 & n_skin == 0 ~ "Nasal-associated",
      TRUE ~ "Other"
    )
  )

top_n_per_habitat <- 15
min_genomes_main <- 5
manual_keep_genera <- c("Dolosigranulum")

skin_ranked_genera <- genus_habitat_counts %>%
  filter(n_skin > 0) %>%
  arrange(desc(n_skin), desc(n_genomes), Genus) %>%
  mutate(
    rank_skin = row_number(),
    in_skin_top = rank_skin <= top_n_per_habitat
  ) %>%
  select(Genus, rank_skin, in_skin_top)

nasal_ranked_genera <- genus_habitat_counts %>%
  filter(n_nasal > 0) %>%
  arrange(desc(n_nasal), desc(n_genomes), Genus) %>%
  mutate(
    rank_nasal = row_number(),
    in_nasal_top = rank_nasal <= top_n_per_habitat
  ) %>%
  select(Genus, rank_nasal, in_nasal_top)

genus_selection_table <- genus_habitat_counts %>%
  left_join(skin_ranked_genera, by = "Genus") %>%
  left_join(nasal_ranked_genera, by = "Genus") %>%
  mutate(
    rank_skin = ifelse(is.na(rank_skin), NA_integer_, rank_skin),
    rank_nasal = ifelse(is.na(rank_nasal), NA_integer_, rank_nasal),
    in_skin_top = ifelse(is.na(in_skin_top), FALSE, in_skin_top),
    in_nasal_top = ifelse(is.na(in_nasal_top), FALSE, in_nasal_top),
    in_ranked_union = in_skin_top | in_nasal_top,
    in_manual_keep = Genus %in% manual_keep_genera,
    include_rebuild_main = (in_ranked_union & n_genomes >= min_genomes_main) | in_manual_keep
  ) %>%
  rowwise() %>%
  mutate(
    selection_reason = paste(
      c(
        if (in_skin_top && n_genomes >= min_genomes_main) paste0("skin_top", top_n_per_habitat) else NULL,
        if (in_nasal_top && n_genomes >= min_genomes_main) paste0("nasal_top", top_n_per_habitat) else NULL,
        if (in_manual_keep) "manual_keep" else NULL
      ),
      collapse = ";"
    ),
    selection_reason = ifelse(selection_reason == "", "not_selected", selection_reason)
  ) %>%
  ungroup()

selected_genera <- genus_selection_table %>%
  filter(include_rebuild_main) %>%
  pull(Genus)

figure3b_selected_genera_metadata <- genus_selection_table %>%
  filter(include_rebuild_main) %>%
  arrange(
    factor(HabitatGroup, levels = c("Skin-associated", "Shared", "Nasal-associated", "Other")),
    desc(n_genomes),
    Genus
  ) %>%
  select(
    Genus,
    HabitatGroup,
    n_skin,
    n_nasal,
    n_genomes,
    n_species,
    rank_skin,
    rank_nasal,
    in_skin_top,
    in_nasal_top,
    in_manual_keep,
    selection_reason
  )

copy_long <- human_surface_wide %>%
  filter(Genus %in% selected_genera) %>%
  select(taxon_oid, Genus, Species, Habitat, all_of(all_cf_ids)) %>%
  pivot_longer(cols = all_of(all_cf_ids), names_to = "CF", values_to = "copy_number")

genome_total_copy <- human_surface_wide %>%
  filter(Genus %in% selected_genera) %>%
  mutate(total_cf_copy = rowSums(across(all_of(all_cf_ids)), na.rm = TRUE)) %>%
  group_by(Genus) %>%
  summarise(mean_total_CF_copy_per_genome = mean(total_cf_copy, na.rm = TRUE), .groups = "drop")

figure3b_genus_cf <- copy_long %>%
  group_by(Genus, CF) %>%
  summarise(
    n_genomes_cf = n(),
    positive_genomes = sum(copy_number > 0, na.rm = TRUE),
    prevalence = positive_genomes / n_genomes_cf,
    mean_copy_all = mean(copy_number, na.rm = TRUE),
    mean_copy_pos = ifelse(positive_genomes > 0, mean(copy_number[copy_number > 0], na.rm = TRUE), 0),
    .groups = "drop"
  ) %>%
  left_join(genus_selection_table, by = "Genus") %>%
  left_join(genome_total_copy, by = "Genus") %>%
  left_join(cf_info %>% select(CF = cf_id, module, functional_theme, functional_theme_short), by = "CF") %>%
  group_by(Genus) %>%
  mutate(CF_richness = sum(prevalence >= 0.05 & positive_genomes >= 2, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(
    main_selection_reason = selection_reason,
    module = as.character(module),
    CF = factor(CF, levels = all_cf_ids),
    module = factor(module, levels = module_order),
    Genus = factor(Genus, levels = genus_selection_table %>%
                     filter(Genus %in% selected_genera) %>%
                     arrange(factor(HabitatGroup, levels = c("Skin-associated", "Shared", "Nasal-associated", "Other")),
                             desc(n_genomes), Genus) %>%
                     pull(Genus))
  ) %>%
  arrange(Genus, module, CF)

species_metrics <- human_surface_wide %>%
  filter(Genus %in% selected_genera) %>%
  group_by(Species, Genus) %>%
  summarise(
    n_genomes = n(),
    n_skin = sum(Habitat == "Skin"),
    n_nasal = sum(Habitat == "Nasal"),
    .groups = "drop"
  )

species_copy_long <- human_surface_wide %>%
  filter(Genus %in% selected_genera) %>%
  select(taxon_oid, Genus, Species, Habitat, all_of(all_cf_ids)) %>%
  pivot_longer(cols = all_of(all_cf_ids), names_to = "CF", values_to = "copy_number")

figure3c_species_cf <- species_copy_long %>%
  group_by(Species, Genus, CF) %>%
  summarise(
    n_genomes_cf = n(),
    positive_genomes = sum(copy_number > 0, na.rm = TRUE),
    prevalence = positive_genomes / n_genomes_cf,
    mean_copy_all = mean(copy_number, na.rm = TRUE),
    mean_copy_pos = ifelse(positive_genomes > 0, mean(copy_number[copy_number > 0], na.rm = TRUE), 0),
    .groups = "drop"
  ) %>%
  left_join(species_metrics, by = c("Species", "Genus")) %>%
  left_join(cf_info %>% select(CF = cf_id, module, functional_theme, functional_theme_short), by = "CF") %>%
  mutate(
    module = factor(as.character(module), levels = module_order),
    CF = factor(CF, levels = all_cf_ids)
  ) %>%
  arrange(Genus, Species, module, CF)

priority_species_table <- tibble::tribble(
  ~Species, ~selection_group, ~selection_reason,
  "Staphylococcus epidermidis", "Skin Staphylococcus", "common epithelial-associated skin Staphylococcus; priority representative",
  "Staphylococcus aureus", "Skin/Nasal Staphylococcus", "clinically important Staphylococcus and nasal carriage species",
  "Staphylococcus hominis", "Skin Staphylococcus", "common epithelial-associated skin Staphylococcus; expanded skin-side representation",
  "Staphylococcus capitis", "Skin Staphylococcus", "common epithelial-associated skin Staphylococcus; expanded skin-side representation",
  "Staphylococcus warneri", "Skin Staphylococcus", "skin-enriched Staphylococcus in current MAG collection",
  "Staphylococcus saprophyticus", "Skin Staphylococcus", "skin-enriched Staphylococcus in current MAG collection",
  "Staphylococcus simulans", "Skin Staphylococcus", "skin-enriched Staphylococcus with n_genomes>=5",
  "Staphylococcus auricularis", "Skin Staphylococcus", "skin-enriched Staphylococcus with n_genomes>=5",
  "Cutibacterium acnes", "Skin Cutibacterium", "canonical epithelial-associated skin Cutibacterium",
  "Cutibacterium granulosum", "Skin Cutibacterium", "epithelial-associated skin Cutibacterium representative",
  "Cutibacterium modestum", "Skin Cutibacterium", "skin-enriched Cutibacterium in current MAG collection",
  "Cutibacterium namnetense", "Skin Cutibacterium", "skin-enriched Cutibacterium with n_genomes>=5",
  "Dolosigranulum pigrum", "Nasal/shared representative", "nasal-relevant commensal representative",
  "Corynebacterium accolens", "Nasal/shared representative", "nasal-associated Corynebacterium representative",
  "Corynebacterium propinquum", "Nasal/shared representative", "nasal-enriched Corynebacterium representative",
  "Corynebacterium pseudodiphtheriticum", "Nasal/shared representative", "nasal-enriched Corynebacterium representative",
  "Moraxella catarrhalis", "Nasal/shared representative", "nasal-associated Moraxella representative",
  "Moraxella nonliquefaciens", "Nasal/shared representative", "nasal-associated Moraxella representative"
)

selected_species_metadata <- priority_species_table %>%
  left_join(species_metrics, by = "Species") %>%
  mutate(
    available_in_current_data = !is.na(n_genomes),
    passes_min_genomes = available_in_current_data & n_genomes >= 5,
    selected_for_figure3c = passes_min_genomes
  )

selected_key_species <- selected_species_metadata %>%
  filter(selected_for_figure3c) %>%
  pull(Species)

figure3c_key_species <- figure3c_species_cf %>%
  filter(Species %in% selected_key_species) %>%
  mutate(Species = factor(Species, levels = selected_key_species)) %>%
  arrange(Species, module, CF)

write.table(figure3a_cf_table, file.path(paths$out_dir, "Figure3A_CF_FunctionalTheme_Table.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(figure3a_module_summary, file.path(paths$out_dir, "Figure3A_Module_FunctionalTheme_Summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(figure3b_selected_genera_metadata, file.path(paths$out_dir, "Figure3B_SelectedGenera_Metadata.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(figure3b_genus_cf, file.path(paths$out_dir, "Figure3B_Genus_AllCF_Deployment.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(figure3c_species_cf, file.path(paths$out_dir, "Figure3C_Species_AllCF_Architecture.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(figure3c_key_species, file.path(paths$out_dir, "Figure3C_KeySpecies_AllCF_Architecture.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(data.frame(SelectedSpecies = selected_key_species), file.path(paths$out_dir, "Figure3C_KeySpecies_Selected.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(selected_species_metadata, file.path(paths$out_dir, "Figure3C_KeySpecies_Rationale.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

rules_md <- c(
  "# Figure3B Genus Selection Rules",
  "",
  "## English",
  "",
  "Genera in Figure 3B were selected as representative genera recurrently recovered from the skin and nasal MAG collection in this study. The selection was data-driven and did not use the previous curated genus-prefix list.",
  "",
  paste0("1. Only MAGs with genus-level annotation and Habitat normalized to Skin or Nasal were considered."),
  paste0("2. Genera were ranked separately by the number of Skin genomes and the number of Nasal genomes."),
  paste0("3. The union of the top ", top_n_per_habitat, " Skin-ranked genera and top ", top_n_per_habitat, " Nasal-ranked genera was retained."),
  paste0("4. Genera selected by ranking were required to contain at least ", min_genomes_main, " genomes in total."),
  paste0("5. Manual keep genera were retained and explicitly flagged: ", paste(manual_keep_genera, collapse = "; "), "."),
  "",
  "This rule supports the wording \"representative genera recovered from the skin and nasal MAG collection\" rather than \"significantly niche-associated genera\".",
  "",
  "## Recommended Caption Wording",
  "",
  "Genera were ranked separately according to the number of genomes recovered from skin and nasal samples. The union of the top-ranked genera from both niches was retained with a minimum total genome count threshold of five genomes. Manually retained genera, if any, were documented in the selection metadata."
)
writeLines(rules_md, file.path(paths$out_dir, "Figure3B_GenusSelection_Rules.md"), useBytes = TRUE)

validation <- data.frame(
  metric = c(
    "n_cf_modules_file",
    "n_modules",
    "n_functional_themes",
    "n_selected_genera",
    "top_n_per_habitat",
    "min_genomes_main",
    "manual_keep_genera",
    "n_genus_cf_rows",
    "expected_genus_cf_rows",
    "n_species",
    "n_key_species",
    "prevalence_min",
    "prevalence_max",
    "contains_Nose_label"
  ),
  value = c(
    length(unique(all_cf_ids)),
    length(unique(as.character(modules_df$module))),
    length(unique(figure3a_cf_table$functional_theme)),
    length(unique(as.character(figure3b_genus_cf$Genus))),
    top_n_per_habitat,
    min_genomes_main,
    paste(manual_keep_genera, collapse = ";"),
    nrow(figure3b_genus_cf),
    length(unique(as.character(figure3b_genus_cf$Genus))) * length(unique(all_cf_ids)),
    length(unique(figure3c_species_cf$Species)),
    length(selected_key_species),
    min(figure3b_genus_cf$prevalence, na.rm = TRUE),
    max(figure3b_genus_cf$prevalence, na.rm = TRUE),
    any(metadata_df$Habitat == "Nose")
  )
)

write.table(validation, file.path(paths$out_dir, "Figure3_rebuild_validation_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

message("Figure3 rebuild data preparation completed.")
message("Output directory: ", paths$out_dir)


