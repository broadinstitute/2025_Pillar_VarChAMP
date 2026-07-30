# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Repository Overview

This repository contains analysis code and figures for the VarChAMP Pillar project publication. It generates manuscript figures from pre-processed data, with external pipelines referenced for full reproducibility.

## Directory Structure

```
2025_Pillar_VarChAMP/
├── 1_allele_collection/           # Allele processing from Zenodo data
├── 2_dms_bms_overlap_analyses/    # DMS/BMS overlap analyses (imaging)
├── 3_integrated_assay_analyses/   # Cross-assay integration
├── utils/                         # Shared utilities (fonts, functions)
├── env.yml                        # Conda environment specification
└── CLAUDE.md                      # This file
```

---

## 1_allele_collection/

Processes and consolidates Pillar allele collection data from Zenodo.

```
1_allele_collection/
├── 1_inputs/
│   └── raw_inputs/
│       └── final_pillar_data_with_clinvar_*.csv.gz   # Raw Zenodo data
├── 2_analyses/
│   └── 0_pillar_allele_collection_process.ipynb      # Main processing notebook
└── 3_outputs/
    └── pillar_snp_alleles.tsv                        # ~92,941 SNP/missense variants
```

**Workflow:**
- Downloads Pillar data from Zenodo (with MD5 validation)
- Filters to SNP/missense variants (single nucleotide changes)
- Standardizes gene names and HGVS protein nomenclature
- Outputs unique variants across 40 genes

---

## 2_dms_bms_overlap_analyses/

Analysis comparing DMS (Deep Mutational Scanning) and BMS (biomarker/phenotype) overlap. Currently contains F9 imaging analysis.

```
2_dms_bms_overlap_analyses/
└── imaging/
    ├── README.md
    ├── 1_inputs/
    │   └── README.md                    # Points to shared inputs
    ├── 2_analyses/
    │   └── F9_analyses/
    │       ├── 0_extract_f9_data.py     # Extract from external pipeline
    │       ├── 1_F9_visualizations.ipynb # Main figures
    │       └── 2_F9_cell_crops.ipynb    # Cell crop images
    ├── 3_outputs/
    │   ├── data/F9/
    │   │   ├── f9_sc_features_minimal.parquet
    │   │   ├── f9_feature_importance.parquet
    │   │   └── f9_variant_zscore_features.parquet
    │   └── pillar_manuscript_figures/
    │       ├── F9_*.svg                 # Main figures
    │       └── cell_crops/*.svg         # Cell images
    └── _external_pipeline/
        └── README.md                    # Links to public snakemake release
```

**External Pipeline:**
- Public release: https://github.com/broadinstitute/2025_varchamp_snakemake/releases/tag/PillarManuscript
- Run batches: 7, 8, 11, 12, 13, 14, 15, 16

**Quick Start:**
```bash
cd 2_dms_bms_overlap_analyses/imaging/2_analyses/F9_analyses/
jupyter notebook 1_F9_visualizations.ipynb
```

---

## 3_integrated_assay_analyses/

Cross-assay integration combining results from DualIPA, PPI, and Imaging assays.

```
3_integrated_assay_analyses/
├── 1_inputs/
│   ├── README.md
│   ├── VarChAMP_data_supp_mat_PP.tsv        # Publication-ready assay results
│   └── 0_all_gene_variants_assayed_summary.tsv  # 1,013-column annotations
├── 2_analyses/
│   └── 0_integrative_assay_summary.ipynb    # Main integration analysis
└── 3_outputs/
    └── pillar_manuscript_figures/
        ├── assay_hit_per_gene.svg
        ├── gene_assayability.svg
        ├── PLP_BLB_hits.svg
        └── assay_phenotyping_score_*.svg
```

**Key Input Files:**

| File | Description |
|------|-------------|
| `VarChAMP_data_supp_mat_PP.tsv` | Authoritative source for assay scores/hits (12 columns) |
| `0_all_gene_variants_assayed_summary.tsv` | Comprehensive annotations (1,013 columns: ClinVar, gnomAD, in silico predictors, G2P DMS data, structural features) |

**Note:** `VarChAMP_data_supp_mat_PP.tsv` is imported from [broadinstitute/2025_laval_submitted](https://github.com/broadinstitute/2025_laval_submitted).

---

## utils/

Shared utility files used across analyses.

```
utils/
├── ARIAL.TTF     # Font for publication figures
└── utils.py      # Shared Python functions
```

**Key Functions in `utils.py`:**
- `plot_assay_hit_by_category()` - Visualization for hit rates by category
- `compute_aubprc()` - Adjusted precision-recall metrics
- `plot_auroc_curves()` - ROC curve analysis
- `plot_gene_level_summary()` - Gene-level hit summaries

---

## Development Environment

### Conda Environment
```bash
source "$HOME/software/anaconda3/etc/profile.d/conda.sh"
conda activate varchamp
```

Or create from spec:
```bash
conda env create -f env.yml
conda activate varchamp
```

### Key Dependencies
- Python 3.8, numpy, pandas, matplotlib, seaborn, scikit-learn
- Polars for high-performance data manipulation
- Jupyter notebook environment

---

## Data Standards

### File Formats
- **TSV/CSV**: Primary data format
- **Parquet**: Intermediate processing (compressed with zstd)
- **Jupyter notebooks**: Analysis workflows
- **SVG**: Publication figures

### Naming Conventions
- Gene alleles: `gene_variant` = `symbol` + "_" + `aa_change` (e.g., `F9_Cys28Arg`)

### Variant classification labels & colors

Figures color variants by a **ClinVar / population-observation** class. These
labels and colors are kept **consistent with the main IGVF paper repo**
([bbi-lab/IGVF-cvfg-pillar-project](https://github.com/bbi-lab/IGVF-cvfg-pillar-project),
Tejura et al. 2026, [bioRxiv 10.64898/2026.02.14.705848](https://doi.org/10.64898/2026.02.14.705848)):

| Label | Definition | Color |
|-------|------------|-------|
| **PLP** | ClinVar Pathogenic / Likely pathogenic | `#CA7682` |
| **BLB** | ClinVar Benign / Likely benign | `#1D7AAB` |
| **Conflicting** | ClinVar conflicting classifications | `#505050` |
| **VUS** | ClinVar "Uncertain significance" **only** | `#A0A0A0` |
| **gnomAD** | Not in ClinVar, observed in gnomAD (has an allele frequency) | `#A0A0A0` |
| **Unobserved** | Absent from **both** ClinVar and gnomAD | `#A0A0A0` |

**IRON RULE — non-ClinVar variants are NEVER labelled "VUS".** `VUS` is reserved
for ClinVar "Uncertain significance". A variant with no ClinVar record is
`gnomAD` (population-observed) or `Unobserved` (absent from ClinVar and gnomAD),
matching the main paper's definitions. The "unknown" family (VUS / gnomAD /
Unobserved) shares the neutral gray `#A0A0A0`; use the lighter neutral `#E0E0E0`
only if two of them must be visually distinguished in the same panel.

- Determine class from `clinvar_clnsig_clean_pp_strict` (ClinVar tier),
  `clinvar_sig_2025` (ClinVar record present?) and `gnomad_MAF` (in gnomAD?).
  Matching must be **case-insensitive** so `Likely pathogenic`/`Likely benign`
  group into PLP/BLB (see the F9 correlation figure in
  `2_dms_bms_overlap_analyses/imaging/2_analyses/F9_analyses/1_F9_visualizations.ipynb`).
- The upstream `clinvar_clnsig_clean*` "Others" bucket is **not** a single
  class: it mixes not-in-ClinVar variants with in-ClinVar-but-non-standard
  calls (`not provided`, `other`, `risk factor`, …). Do not equate "Others"
  with "Unobserved" — resolve using ClinVar-record + gnomAD presence.

---

## Running Analyses

### Imaging Analysis (F9)
```bash
cd 2_dms_bms_overlap_analyses/imaging/2_analyses/F9_analyses/

# Visualizations use pre-extracted parquet files
jupyter notebook 1_F9_visualizations.ipynb
jupyter notebook 2_F9_cell_crops.ipynb

# To regenerate parquet files from external pipeline outputs:
python 0_extract_f9_data.py
```

### Integrated Analysis
```bash
cd 3_integrated_assay_analyses/2_analyses/
jupyter notebook 0_integrative_assay_summary.ipynb
```

### Allele Collection
```bash
cd 1_allele_collection/2_analyses/
jupyter notebook 0_pillar_allele_collection_process.ipynb
```

---

## External Resources

| Resource | URL |
|----------|-----|
| VarChAMP Snakemake Pipeline | https://github.com/broadinstitute/2025_varchamp_snakemake/releases/tag/PillarManuscript |
| Laval Submitted (assay results) | https://github.com/broadinstitute/2025_laval_submitted (still private, will be released soon with the paper) |
