# Microbial Colonization Modules in Skin and Nasal Microbiomes

This repository contains the source code for plotting and analyzing the genetic basis of microbial colonization in human skin and nasal microbiomes, as described in our study.

## Project Overview

Microbial colonization is essential for long-term persistence in host-associated ecosystems. This project systematically investigates the genetic organization of colonization determinants in exposed epithelial habitats (skin and nasal) using a large-scale comparative genomics framework.

### Key Findings
- **High-Quality Resource**: Constructed a reference collection of 7,248 metagenome-assembled genomes (MAGs) across nasal, skin, and environmental habitats.
- **Colonization Factors (CFs)**: Identified 31 CFs through Genotype–Habitat Association (GHA) analysis, focusing on the top 0.2% of protein families exhibiting strongest host enrichment.
- **Colonization Modules (CMs)**: Resolved 18 CMs reflecting coordinated genetic systems for nutrient acquisition (e.g., methionine uptake), environmental stress adaptation (e.g., oxidative stress tolerance), and cellular maintenance.
- **Organizational Principles**: Demonstrated that ecological specialization, rather than phylogenetic relatedness, drives the diversification and deployment of colonization systems on epithelial surfaces.

## Repository Structure

The repository is organized to facilitate the reproduction of figures and analyses presented in the manuscript.

```text
paper_submission_code/
├── scripts/
│   ├── plotting/             # R scripts for generating manuscript figures
│   │   ├── Figure1/          # Dataset construction, quality control, and taxonomy
│   │   ├── Figure2/          # GHA analysis, CF identification, and modularity
│   │   ├── Figure3/          # Taxonomic deployment and functional landscapes
│   │   └── Figure4/          # Evolutionary distribution and GTDB-scale expansion
│   └── utils/
│       └── 配色.R             # Universal color configuration and themes
├── data/                     # [Placeholder] Processed data for plotting
└── README.md                 # Project documentation
```

## Getting Started

### Prerequisites
To run the plotting scripts, you will need **R (>= 4.0.0)** with the following packages:
- `tidyverse`, `ggplot2`, `patchwork`, `scales`, `vroom`, `RColorBrewer`, `ComplexHeatmap`

### Usage
Each figure directory contains the specific R scripts required to generate the corresponding panels. For example, to generate Figure 1:
```R
source("scripts/plotting/Figure1/Figure1C_Taxonomic_Diversity.R")
```
*Note: Ensure the working directory is set to the root of this repository. Most scripts depend on `scripts/utils/配色.R` for consistent visual styling.*

## Data Availability
The genomic data used in this study include MAGs from the 4D-SZ cohort, UHSG, HSMG, and the GEM catalogue. Processed phylogenetic profiles and MI Z-scores are available in the `data/` directory (or via the link provided in the manuscript).

## Citation
If you use this code or our findings in your research, please cite:
> [Author Names]. (2026). Systems-level identification of microbial colonization modules in human skin and nasal microbiomes. *[Journal Name]*.

---
For questions or issues, please open an issue on this GitHub repository.
