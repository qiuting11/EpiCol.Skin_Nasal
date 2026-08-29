# Data availability and external file manifest

This repository includes small derived tables under `data/figure1` through `data/figure5`. Large raw or intermediate files are intentionally excluded from the GitHub draft.

To run scripts that require excluded inputs, place the files under `data/external/` with the names shown below, or edit the corresponding script paths.

## Excluded large inputs

| Expected path in repository | Original local source | Approx. size | Used by | Reason excluded |
| --- | --- | ---: | --- | --- |
| `data/external/mi_zscore.csv` | `E:/data/gephe_output_11/association/result_11_1230/mi_zscore.csv` | 1.2 GB | `Figure2A_MI_Zscore_Distribution.R` | Too large for GitHub; full MI Z-score table. |
| `data/external/eggnog_output_with_pog_multi_level.tsv` | `E:/data/gephe_output_11/analysis/input/eggnog_output_with_pog_multi_level.tsv` | 12.9 MB | `Figure3_rebuild_prepare_data.R` | Optional annotation rebuild input; omitted to keep repository compact. |
| `data/external/figure2_homology_network/pog_relationships.tsv` and `data/external/figure2_homology_network/pog_info.tsv` | upstream POG homology-network working tables | intermediate | `Figure2C_CF_Homology_Network.R` | Homology-network intermediate tables; omitted from lightweight GitHub release. |
| `data/external/figure3D_three_genus_pcoa/mag_cf_long_with_core_phenotypes.tsv` | `E:/data/gephe_output_11/human_hq_mag_cf_matrix/main_strict_sensitive_result/visualization_analysis/tables/mag_cf_long_with_core_phenotypes.tsv` | intermediate | `Figure3D_three_genus_pcoa_vegan.R` | MAG-level eCF copy-number table used to build all-HQ-MAG three-genus PCoA matrices. |
| `data/external/true_cf_genome_copy_counts.csv` | `E:/data/mmseqs/sp/04_Final_Metrics/true_cf_genome_copy_counts.csv` | 37.2 MB | `Figure4C_Dosage_Expansion.R` | GTDB-scale copy-count table; larger derived matrix. |
| `data/external/bac120_taxonomy.tsv` | `E:/data/mmseqs/sp/02_Reference/bac120_taxonomy.tsv` | 102.6 MB | `Figure4C_Dosage_Expansion.R` | GTDB taxonomy table; external reference-derived data. |
| `data/external/supplementary_s2_s3/top0.2_I3.0_modules.txt` | `E:/data/gephe_output_11/analysis/output_test_1/top0.2_I3.0_modules.txt` | small | `Supplementary_S2_S3/Fig_S2D_to_S2G_conservative_filtering_delta.py` | eCF/eCM order for current supplementary S2/S3 panels. |
| `data/external/supplementary_s2_s3/Figure3_Final/*.tsv` | `E:/data/R/gephe_R/results/figures/Figure3_Final/*.tsv` | small/intermediate | `Supplementary_S2_S3/Fig_S2D_to_S2G_conservative_filtering_delta.py` | Previous broad homolog-recovery display tables used for conservative-filtering deltas. |
| `data/external/supplementary_s2_s3/Figure3_human_pre11_paper_id50_k1_original_style/*.tsv` | `E:/data/R/gephe_R/results/figures/Figure3_human_pre11_paper_id50_k1_original_style/*.tsv` | small/intermediate | `Supplementary_S2_S3/Fig_S2B*`, `Fig_S2C*`, `Fig_S2D_to_S2G*` | Epithelial-associated MAG subset display tables and selected taxon lists. |
| `data/external/supplementary_s2_s3/human_hq_mag_cf_matrix/controlled_parameter_comparison/tables/*.tsv` | `E:/data/gephe_output_11/human_hq_mag_cf_matrix/controlled_parameter_comparison/tables/*.tsv` | small/intermediate | `Supplementary_S2_S3/Fig_S3A*` | Controlled DIAMOND parameter comparison summary tables. |
| `data/external/supplementary_s2_s3/human_hq_mag_cf_matrix/main_strict_sensitive_result/cf_presence_absence_analysis/tables/*.tsv` | `E:/data/gephe_output_11/human_hq_mag_cf_matrix/main_strict_sensitive_result/cf_presence_absence_analysis/tables/*.tsv` | small/intermediate | `Supplementary_S2_S3/Fig_S3A*`, `Fig_S3B*`, `Fig_S3C_to_S3F*` | Strict-sensitive presence/absence association and plotting tables. |
| `data/external/supplementary_s2_s3/human_hq_mag_cf_matrix/sensitivity_*/outputs/mag_cf_presence_with_phenotypes.tsv` | `E:/data/gephe_output_11/human_hq_mag_cf_matrix/sensitivity_*/outputs/mag_cf_presence_with_phenotypes.tsv` | intermediate | `Supplementary_S2_S3/Fig_S3B*` | MAG-paired presence matrices for McNemar tests across controlled parameter settings. |

## Included small derived data

- Figure 1: phenotype metadata, selected genome metadata, CheckM2 summaries, and source data tables needed by Figure 1 and Figure S1 plotting scripts.
- Figure 2: MCL cluster summaries, top eCF/eCM module tables, selected co-cluster data, and compact phylogenetic profile input.
- Figure 3: module metadata, selected genus/species deployment tables, copy-number profiles, and grand butterfly-map input.
- Figure 4: reduced GTDB phylum tree, phylum sizes, eCF phylum coverage, and module table.
- Figure 5: compact synteny tables for metNIQ and serine dehydratase architecture plots.
- Supplementary S2/S3: script-level manifest and final panel rename manifest are included; working-analysis input tables are treated as external inputs.

## Reproducibility note

Some scripts can run fully using the bundled `data/` inputs. Scripts that rebuild analyses from large upstream matrices require the excluded files above. The repository is therefore a clean code-and-small-data draft, not a complete raw-data archive.
