# Epithelial-associated microbial eCF/eCM analysis code

This repository contains the cleaned plotting code and small derived data tables used to reproduce the manuscript figure panels for an epithelial-associated versus environmental microbial genome analysis.

The repository intentionally includes only code and small processed/derived tables. Large raw matrices, full genome resources, and external catalogues are not bundled; see `docs/data_availability.md`.

## Repository layout

```text
paper_submission_code/
  data/                       # small processed data included in this repository
    figure1/                  # Figure 1 and Figure S1 plotting inputs
    figure2/                  # Figure 2 small plotting inputs
    figure3/                  # Figure 3 small plotting inputs and intermediate tables
    figure4/                  # Figure 4 small plotting inputs
    figure5/                  # Figure 5 synteny plotting inputs and eCM1/eCM13 AF input files
  docs/
    data_availability.md      # external/large data manifest
    repository_manifest.tsv   # tracked-file inventory with sizes
    tree.txt                  # tracked directory tree
  outputs/                    # generated figures; ignored by git except .gitkeep
    Figure1/
    Figure2/
    Figure3/
    Figure4/
    Figure5/
    Supplementary_S2_S3/
  scripts/
    plotting/
      Figure1/
      Figure2/
      Figure3/
      Figure4/
      Figure5/
      Supplementary_S2_S3/
    utils/
  README.md
```

## Terminology used in figures

- `eCF`: ecological Candidate Factor. Internal data columns may still use `CF`, `CF_id`, or `cf_id` for compatibility.
- `eCM`: ecological Candidate Module. Internal data columns may still use `module`.
- `Epithelial-associated`: used for microbial genomes recovered from epithelial habitats such as skin microenvironments and nasal cavity.

## Requirements

Tested with R 4.5.x. Required R packages include:

```r
install.packages(c(
  "tidyverse", "ggplot2", "dplyr", "tidyr", "readr", "vroom",
  "patchwork", "scales", "stringr", "forcats", "ggrepel",
  "RColorBrewer", "igraph", "circlize", "ComplexHeatmap", "ape",
  "ggnewscale"
))
```

Some Figure 4 scripts also use Bioconductor/YuLab packages:

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("ggtree", "ggtreeExtra"))
```

## How to run

Run scripts from the repository root:

```r
setwd("path/to/paper_submission_code")
source("scripts/plotting/Figure1/Figure1B_Selection_Efficiency_Modern.R")
```

Generated figures are written to `outputs/Figure1`, `outputs/Figure2`, `outputs/Figure3`, `outputs/Figure4`, or `outputs/Figure5`.

Recommended order for Figure 3:

```r
source("scripts/plotting/Figure3/Figure3_rebuild_prepare_data.R")
source("scripts/plotting/Figure3/Figure3B_Genus_Deployment.R")
source("scripts/plotting/Figure3/FigureS2A_Species_Architecture.R")
source("scripts/plotting/Figure3/FigureS2B_Grand_Butterfly_Map.R")
source("scripts/plotting/Figure3/Figure3D_three_genus_pcoa/Figure3D_three_genus_pcoa_vegan.R")
```


Figure 5 synteny panels:

```r
source("scripts/plotting/Figure5/Figure5_metNIQ_real_coords.R")
source("scripts/plotting/Figure5/Figure5_sda_synteny.R")
```
The eCM1 and eCM13 AlphaFold 3 input files are stored under `data/figure5/eCM1_AF_input` and `data/figure5/eCM13_AF_input`, respectively.


Current supplementary Figure S2/S3 panels:

```r
source("scripts/plotting/Supplementary_S2_S3/run_supplementary_s2_s3.R")
```

These scripts require excluded working-analysis inputs under `data/external/supplementary_s2_s3/`; see `scripts/plotting/Supplementary_S2_S3/README.md` and `docs/data_availability.md`. Figure-facing labels use `Strict-fast`, `Strict-sensitive`, `Broad homolog-recovery`, `Cutibacterium`, `eCF`, and `eCM`. Legacy/internal data values are only retained as mapping keys where needed to read existing tables.

Scripts that require large excluded data files will stop with a missing-file error unless those files are supplied under data/external/; see docs/data_availability.md.

The Fig. 3D three-genus PCoA script expects `data/external/figure3D_three_genus_pcoa/mag_cf_long_with_core_phenotypes.tsv` unless an explicit input path is supplied.

## Included data policy

Included files are small derived tables needed for plotting and review. Excluded files include large MI Z-score tables, GTDB-wide copy-count matrices, full GTDB taxonomy, homology-network intermediates, supplementary working-analysis inputs, generated figures, and genome sequences. These are documented in `docs/data_availability.md`.

## Validation performed for this draft

- All `.R` files under `scripts/` parse successfully with R.
- Manuscript-facing scripts use repository-relative paths; external inputs are documented under `docs/data_availability.md`.
- FigureS2B filters out the `CF22/CF25` panel (`Pattern 5 (RNR Alpha)`) and keeps the remaining four panels.

## Citation

If you use this repository, please cite the associated manuscript once available.


