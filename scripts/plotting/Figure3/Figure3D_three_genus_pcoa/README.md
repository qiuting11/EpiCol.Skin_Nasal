# Figure 3D three-genus eCF PCoA

This directory contains scripts for the Fig. 3D candidate PCoA analysis of eCF profiles across all HQ MAGs from three representative epithelial-associated genera.

## Scope

Included genera:

- `Staphylococcus`
- `Corynebacterium`
- `Cutibacterium_or_Propionibacterium`, displayed as `Cutibacterium`

This is an all-HQ-MAG three-genus analysis. It is not the `human_pre11_paper_id50_k1` conservative subset.

## Primary script

Use the R script for manuscript-facing analysis and plotting:

```r
source("scripts/plotting/Figure3/Figure3D_three_genus_pcoa/Figure3D_three_genus_pcoa_vegan.R")
```

Default input path:

```text
data/external/figure3D_three_genus_pcoa/mag_cf_long_with_core_phenotypes.tsv
```

Default output path:

```text
outputs/Figure3/Figure3D_three_genus_pcoa/
```

You can also run with explicit paths:

```bash
Rscript scripts/plotting/Figure3/Figure3D_three_genus_pcoa/Figure3D_three_genus_pcoa_vegan.R path/to/mag_cf_long_with_core_phenotypes.tsv outputs/Figure3/Figure3D_three_genus_pcoa
```

The number of permutations defaults to 999 and can be changed with `N_PERM`.

## Method

Presence/absence profile:

- `presence_binary = copy_number > 0`
- `vegan::vegdist(method = "jaccard", binary = TRUE)`
- `vegan::adonis2()` for PERMANOVA
- `vegan::betadisper()` and `permutest()` for dispersion testing

Copy-number profile:

- `log1p(copy_number)`
- `vegan::vegdist(method = "bray")`
- `vegan::adonis2()` for PERMANOVA
- `vegan::betadisper()` and `permutest()` for dispersion testing

The script also writes negative eigenvalue diagnostics from the principal-coordinate representation used by `betadisper`.

## Outputs

The R script writes:

- `Figure3D_three_genus_presence_jaccard_pcoa.pdf/png`
- `Figure3D_three_genus_presence_jaccard_pcoa_coordinates.tsv`
- `Figure3D_three_genus_presence_jaccard_pcoa_negative_eigen_summary.tsv`
- `Figure3D_three_genus_copy_braycurtis_pcoa.pdf/png`
- `Figure3D_three_genus_copy_braycurtis_pcoa_coordinates.tsv`
- `Figure3D_three_genus_copy_braycurtis_pcoa_negative_eigen_summary.tsv`
- `Figure3D_three_genus_pcoa_statistics.tsv`
- `Figure3D_three_genus_pcoa_sample_counts.tsv`


## Manuscript text

Recommended Nature-family caption and Methods wording is provided in:

`docs/Figure3D_three_genus_pcoa_caption_methods.md`

## Legacy exploratory script

`run_three_genus_ecf_pcoa_python_exploratory.py` is the earlier Python implementation used for exploratory plotting and cross-checking. The R `vegan` script should be treated as the primary manuscript-facing workflow.
