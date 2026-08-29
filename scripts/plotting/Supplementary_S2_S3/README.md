# Supplementary Figures S2-S3 plotting scripts

This directory contains the submission-facing scripts for the current supplementary Figure S2/S3 set. The scripts preserve the final figure logic used in the working analysis, but are organized to run from the `paper_submission_code` repository root.

## Scope

- Fig. S2B-S2G: epithelial-associated MAG eCF/eCM deployment and conservative-filtering delta panels.
- Fig. S3A: controlled DIAMOND parameter comparison bubble plot by genus.
- Fig. S3B: paired McNemar delta-prevalence heatmap by genus and eCF.
- Fig. S3C-S3F: strict-sensitive presence/absence taxonomic association visualizations.
- Fig. S4 is intentionally not included in this round.

## External Input Layout

Large/intermediate inputs are not bundled in the submission-code draft. To rerun these scripts, place the required working-analysis files under:

```text
data/external/supplementary_s2_s3/
  top0.2_I3.0_modules.txt
  Figure3_Final/
    Figure3B_Genus_AllCF_Deployment.tsv
    Figure3C_KeySpecies_AllCF_Architecture_v4.tsv
  Figure3_human_pre11_paper_id50_k1_original_style/
    Figure3B_Genus_AllCF_Deployment.tsv
    Figure3B_SelectedGenera_Metadata.tsv
    Figure3C_Species_AllCF_Architecture.tsv
    Figure3C_KeySpecies_AllCF_Architecture_v4.tsv
    Figure3C_KeySpecies_Selected_v4.tsv
  human_hq_mag_cf_matrix/
    controlled_parameter_comparison/tables/*.tsv
    main_strict_sensitive_result/cf_presence_absence_analysis/tables/*.tsv
    sensitivity_strict_fast/outputs/mag_cf_presence_with_phenotypes.tsv
    sensitivity_strict_sensitive/outputs/mag_cf_presence_with_phenotypes.tsv
    sensitivity_paper_like/outputs/mag_cf_presence_with_phenotypes.tsv
```

The original local sources are documented in `docs/data_availability.md`.

## Run Order

Run from repository root:

```r
source("scripts/plotting/Supplementary_S2_S3/Fig_S2B_genus_eCF_eCM_deployment.R")
source("scripts/plotting/Supplementary_S2_S3/Fig_S2C_species_eCF_eCM_architecture.R")
source("scripts/plotting/Supplementary_S2_S3/Fig_S3A_controlled_parameter_eCF_bubble.R")
source("scripts/plotting/Supplementary_S2_S3/Fig_S3B_controlled_parameter_mcnemar_delta.R")
source("scripts/plotting/Supplementary_S2_S3/Fig_S3C_to_S3F_strict_sensitive_presence_taxonomy.R")
```

Or run the wrapper:

```r
source("scripts/plotting/Supplementary_S2_S3/run_supplementary_s2_s3.R")
```

Outputs are written to `outputs/Supplementary_S2_S3/`.

## Statistical Notes

Presence is defined as `copy_number > 0`. In bubble plots, point size is eCF-positive MAG prevalence, while color shows mean copy number among eCF-positive MAGs as auxiliary information.

For Fig. S3B, paired McNemar tests compare the same MAGs between parameter settings for each genus-eCF cell. FDR correction is applied separately within each contrast. Asterisks mark McNemar FDR < 0.05.
