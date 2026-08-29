# eCM13 AF input

This directory contains AlphaFold 3 input files prepared for eCM13/SDA targets shown in Figure 5.

- `generate_ecm13_af3_json.py`: extracts the Figure 5 eCM13 SDA targets and writes AF3 inputs.
- `eCM13_sda_sequences.fasta`: nine target protein sequences extracted from `m13_targets_combined.faa`.
- `eCM13_sda_AF3_batch.json`: AF3 batch input JSON with one job per MAG.

Split SDA architectures are represented as two-chain jobs (`sdaAB` beta + `sdaAA` alpha). Fused SDA architectures are represented as one-chain jobs (`sdaA`). Terminal `*` characters are removed from protein sequences before writing AF3 inputs.

To regenerate from the repository root, provide the external source FASTA:

```bash
python data/figure5/eCM13_AF_input/generate_ecm13_af3_json.py --source-fasta path/to/m13_targets_combined.faa
```
