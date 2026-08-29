# Figure 3D Three-Genus PCoA Caption and Methods Text

Purpose: manuscript-facing caption and Methods language for a Nature-family submission. This text corresponds to the all-HQ-MAG three-genus eCF profile PCoA, not the `human_pre11_paper_id50_k1` conservative subset.

## Recommended Figure Caption Text

### If Fig. 3D shows only the presence/absence PCoA

**(D) Genus-associated structure of genome-level eCF deployment profiles across representative epithelial-associated lineages.** Principal coordinate analysis (PCoA) was performed using Jaccard distances calculated from binary eCF presence-absence profiles across all high-quality MAGs assigned to *Staphylococcus*, *Corynebacterium* and *Cutibacterium*. Each point represents one MAG and colors denote genus. Statistical significance was assessed by PERMANOVA using `vegan::adonis2` with 999 permutations, and differences in within-genus dispersion were evaluated using `vegan::betadisper`. Genus explained a large fraction of eCF profile dissimilarity (R2 = 0.832, P <= 0.001), although dispersion also differed significantly among genera (betadisper, P <= 0.001), indicating that both genus-level separation and within-genus heterogeneity contributed to the ordination structure.

### If Fig. 3D shows both presence/absence and copy-number PCoA panels

**(D) Genus-associated structure of genome-level eCF deployment and copy-number profiles across representative epithelial-associated lineages.** Principal coordinate analysis (PCoA) was performed for all high-quality MAGs assigned to *Staphylococcus*, *Corynebacterium* and *Cutibacterium*. Left, Jaccard distances calculated from binary eCF presence-absence profiles. Right, Bray-Curtis distances calculated from log1p-transformed eCF copy-number profiles. Each point represents one MAG and colors denote genus. Statistical significance was assessed by PERMANOVA using `vegan::adonis2` with 999 permutations, and differences in within-genus dispersion were evaluated using `vegan::betadisper`. Genus explained a large fraction of profile dissimilarity for both presence-absence profiles (R2 = 0.832, P <= 0.001) and copy-number profiles (R2 = 0.789, P <= 0.001). Significant dispersion differences were also observed in both analyses (betadisper, P <= 0.001), indicating that within-genus heterogeneity contributed to the observed ordination structure.

## Recommended Main-Text Sentence

To further evaluate whether eCF repertoires vary coherently across epithelial-associated lineages, we ordinated genome-level eCF profiles from all high-quality MAGs assigned to three representative genera. Both binary deployment profiles and log-transformed copy-number profiles showed strong genus-associated structure, with genus explaining 83.2% and 78.9% of the respective distance variation by PERMANOVA, although significant dispersion differences indicated additional within-genus heterogeneity.

## Recommended Methods Text

### Principal coordinate analysis of eCF profiles across representative genera

To assess whether epithelial colonization factor (eCF) repertoires exhibit genus-associated deployment patterns, we analysed all high-quality MAGs assigned to three representative epithelial-associated genera: *Staphylococcus*, *Corynebacterium* and *Cutibacterium*. MAG-level eCF copy numbers were obtained from the strict-sensitive human HQ MAG eCF matrix. The internal group label `Cutibacterium_or_Propionibacterium` was displayed as *Cutibacterium* for figure presentation. The analysis included 2,053 MAGs in total: 501 *Staphylococcus*, 532 *Corynebacterium* and 1,020 *Cutibacterium* MAGs.

Two complementary eCF profile representations were evaluated. For the deployment analysis, eCF copy numbers were converted to binary presence-absence profiles, with an eCF considered present in a MAG when its copy number was greater than zero. Pairwise dissimilarities among MAGs were then calculated using Jaccard distance. For the copy-number analysis, raw eCF copy numbers were log1p-transformed and pairwise dissimilarities were calculated using Bray-Curtis distance. Principal coordinate analysis (PCoA) was used to visualize the resulting distance matrices.

Genus-associated differences in eCF profile composition were tested by PERMANOVA using `adonis2` in the R package `vegan`, with 999 permutations. Because PERMANOVA can be sensitive to differences in within-group dispersion, homogeneity of multivariate dispersion was evaluated using `betadisper` followed by permutation testing with 999 permutations. Negative eigenvalue diagnostics from the principal-coordinate representation were recorded for both distance matrices. P values from permutation tests are reported as P <= 0.001 when no permuted statistic was as extreme as the observed statistic, corresponding to the minimum attainable P value with 999 permutations.

## Reporting Notes

- Use the presence-absence/Jaccard PCoA as the preferred Fig. 3D panel if the main claim is eCF deployment.
- Use the copy-number/Bray-Curtis PCoA as an additional panel or supplementary sensitivity analysis if space allows.
- Report `P <= 0.001`, not simply `P = 0.001`, because 999 permutations impose a lower bound of 0.001.
- Do not state that PCoA itself proves significance; PCoA is the visualization, whereas PERMANOVA provides the statistical test.
- Because betadisper is significant, avoid wording that attributes the full pattern only to centroid separation among genera.
- A cautious phrase is: "genus-associated structure, with significant within-genus dispersion".

## R-Based Results Used For Text

| Profile | Distance | MAGs | eCFs | PCoA1 | PCoA2 | PERMANOVA R2 | PERMANOVA P | betadisper P | Negative eigenvalue fraction |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| Presence-absence | Jaccard | 2053 | 31 | 44.9% | 28.1% | 0.832 | <=0.001 | <=0.001 | 0.158 |
| Copy number | Bray-Curtis after log1p | 2053 | 31 | 36.9% | 24.7% | 0.789 | <=0.001 | <=0.001 | 0.216 |