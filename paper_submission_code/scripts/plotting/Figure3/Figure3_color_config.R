# Shared color and annotation configuration for Figure 3 rebuild.

full_module_colors <- c(
  "module1"  = "#72190E",
  "module2"  = "#f94144",
  "module3"  = "#F7557F",
  "module4"  = "#7C417F",
  "module5"  = "#A05992",
  "module6"  = "#f3722c",
  "module7"  = "#f8961e",
  "module8"  = "#f9c74f",
  "module9"  = "#F7F056",
  "module10" = "#C2B923",
  "module11" = "#A9D88C",
  "module12" = "#90be6d",
  "module13" = "#43aa8b",
  "module14" = "#277da1",
  "module15" = "#1BA3C6",
  "module16" = "#B9F2F0",
  "module17" = "#4F7CBA",
  "module18" = "#5289C7"
)

functional_theme_map <- c(
  "module1"  = "Nutrient acquisition and substrate utilization",
  "module3"  = "Nutrient acquisition and substrate utilization",
  "module5"  = "Nutrient acquisition and substrate utilization",
  "module9"  = "Nutrient acquisition and substrate utilization",
  "module13" = "Nutrient acquisition and substrate utilization",
  "module2"  = "Nucleotide metabolism and genome maintenance",
  "module6"  = "Nucleotide metabolism and genome maintenance",
  "module8"  = "Nucleotide metabolism and genome maintenance",
  "module4"  = "Stress adaptation and envelope homeostasis",
  "module14" = "Stress adaptation and envelope homeostasis",
  "module10" = "Translation and RNA-centered regulation",
  "module11" = "Translation and RNA-centered regulation",
  "module12" = "Translation and RNA-centered regulation",
  "module15" = "Translation and RNA-centered regulation",
  "module18" = "Intercellular communication signal",
  "module7"  = "Poorly characterized candidate modules",
  "module16" = "Poorly characterized candidate modules",
  "module17" = "Poorly characterized candidate modules"
)

functional_theme_colors <- c(
  "Nutrient acquisition and substrate utilization" = "#4C78A8",
  "Nucleotide metabolism and genome maintenance" = "#72B7B2",
  "Stress adaptation and envelope homeostasis" = "#E45756",
  "Translation and RNA-centered regulation" = "#54A24B",
  "Intercellular communication signal" = "#B279A2",
  "Poorly characterized candidate modules" = "#9D755D"
)

functional_theme_short <- c(
  "Nutrient acquisition and substrate utilization" = "Nutrient/substrate use",
  "Nucleotide metabolism and genome maintenance" = "Nucleotide/genome maintenance",
  "Stress adaptation and envelope homeostasis" = "Stress/envelope homeostasis",
  "Translation and RNA-centered regulation" = "Translation/RNA regulation",
  "Intercellular communication signal" = "Communication signal",
  "Poorly characterized candidate modules" = "Poorly characterized"
)

functional_theme_relevance <- c(
  "Nutrient acquisition and substrate utilization" = "Resource acquisition under nutrient-limited surface conditions",
  "Nucleotide metabolism and genome maintenance" = "Genome maintenance and salvage under constrained growth",
  "Stress adaptation and envelope homeostasis" = "Response to oxidative, osmotic, and envelope stresses",
  "Translation and RNA-centered regulation" = "Growth tuning and post-transcriptional regulation",
  "Intercellular communication signal" = "Potential community-level signaling",
  "Poorly characterized candidate modules" = "Candidate mechanisms requiring further validation"
)

module_order <- paste0("module", seq_len(18))
