# Imaging Analysis

This directory contains the Cell Painting imaging analysis for the Pillar project, focusing on F9 variant visualization.

## Directory Structure

```
imaging/
├── README.md                          # This file
├── 1_inputs/
│   └── README.md                      # Points to shared inputs
├── 2_analyses/
│   └── F9_analyses/                   # Active F9 analysis
│       ├── 0_extract_f9_data.py       # Extract data from snakemake outputs
│       ├── 1_F9_visualizations.ipynb  # Main figure generation
│       └── 2_F9_cell_crops.ipynb      # Cell crop images
├── 3_outputs/
│   ├── pillar_manuscript_figures/     # All SVG figures
│   │   └── cell_crops/
│   └── data/F9/                       # Extracted parquet files
└── _external_pipeline/
    └── README.md                      # Links to public snakemake release
```

## Quick Start

The visualization notebooks use pre-extracted parquet files (~55MB) that are already in the repository:

```bash
cd 2_analyses/F9_analyses/
jupyter notebook 1_F9_visualizations.ipynb
jupyter notebook 2_F9_cell_crops.ipynb
```

## Data Flow

```
[External Snakemake Pipeline]          See _external_pipeline/README.md
         │
         v
3_integrated_assay_analyses/1_inputs/  Shared input annotations
         │
         v
0_extract_f9_data.py                   Extract minimal F9 data
         │
         v
3_outputs/data/F9/*.parquet            Pre-extracted files (~55MB)
         │
         v
1_F9_visualizations.ipynb              Generate manuscript figures
2_F9_cell_crops.ipynb                  Generate cell crop images
         │
         v
3_outputs/pillar_manuscript_figures/   Final SVG figures
```

## Reproducing from Scratch

To regenerate the parquet files from raw pipeline outputs:

1. Run the [VarChAMP snakemake pipeline](https://github.com/broadinstitute/2025_varchamp_snakemake/releases/tag/PillarManuscript) for batches 7, 8, 11, 12, 13, 14, 15, 16

2. Update paths in `0_extract_f9_data.py` to point to your outputs

3. Run extraction:
   ```bash
   cd 2_analyses/F9_analyses/
   python 0_extract_f9_data.py
   ```

## Output Figures

All figures are saved to `3_outputs/pillar_manuscript_figures/`:
- `F9_MAVE_Img_Corr.svg` - Imaging vs MAVE correlation
- `F9_sc_gfp_feat_dist_main.svg` - Main violin plots
- `F9_sc_gfp_feat_dist_supp.svg` - Supplementary violin plots
- `F9_variant_heatmap.svg` - Z-score feature heatmap
- `cell_crops/*.svg` - Individual cell crop images

## Variant classification labels (ClinVar / observation)

The `F9_MAVE_Img_Corr.svg` scatter colors variants by ClinVar / population
class, kept **consistent with the main IGVF paper repo**
([bbi-lab/IGVF-cvfg-pillar-project](https://github.com/bbi-lab/IGVF-cvfg-pillar-project),
Tejura et al. 2026):

| Label | Definition | Color |
|-------|------------|-------|
| PLP | ClinVar Pathogenic / Likely pathogenic | `#CA7682` |
| BLB | ClinVar Benign / Likely benign | `#1D7AAB` |
| Conflicting | ClinVar conflicting classifications | `#505050` |
| VUS | ClinVar "Uncertain significance" only | `#A0A0A0` |
| gnomAD | Not in ClinVar, observed in gnomAD | `#A0A0A0` |
| Unobserved | Absent from both ClinVar and gnomAD | `#A0A0A0` |

**Non-ClinVar variants are never labelled "VUS."** For F9 the two variants with
no ClinVar record (`Cys407Arg`, `Glu73Lys`) are also absent from gnomAD, so they
are shown as **Unobserved** (not VUS). See the full convention and rationale in
the repo-root [`CLAUDE.md`](../../CLAUDE.md) → *Data Standards → Variant
classification labels & colors*.

## Environment

Requires the `varchamp` conda environment:
```bash
source "$HOME/software/anaconda3/etc/profile.d/conda.sh"
conda activate varchamp
```
