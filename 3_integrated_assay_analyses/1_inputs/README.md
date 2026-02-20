## VarChAMP_data_supp_mat_PP.tsv

The input file [VarChAMP_data_supp_mat_PP.tsv](https://github.com/broadinstitute/2025_laval_submitted/blob/main/3_integrated_assay_results/3_outputs/VarChAMP_data_supp_mat_PP.tsv) is imported from https://github.com/broadinstitute/2025_laval_submitted with commit @532d860a938eacc86055158b099b99352b5d3339.

## 0_all_gene_variants_assayed_summary.tsv

This file contains comprehensive annotations for all 1,407 variants assayed in VarChAMP. The variants in this file are aligned with `VarChAMP_data_supp_mat_PP.tsv` (matched by `symbol` + `aa_change`). The file contains 1,013 columns organized into the following categories:

**Note:** VarChAMP assay results (DualIPA, PPI, Imaging scores and hits) are NOT included in this file. For assay results, use `VarChAMP_data_supp_mat_PP.tsv` which contains the authoritative publication-ready values.

### Variant Identification (Columns 1-16)
| Column | Description |
|--------|-------------|
| `symbol` | Gene symbol (e.g., BRCA1, TP53) |
| `aa_change` | Amino acid change (e.g., Arg389Gly) |
| `gene_variant` | Combined gene_aa_change identifier |
| `spdi` | SPDI notation for the variant |
| `nt_change` | Nucleotide change |
| `chr_num` | Chromosome number |
| `nuc_loc` | Nucleotide location |
| `ref_allele` | Reference allele |
| `alt_allele` | Alternate allele |
| `Chrom` | Chromosome |
| `ensembl_gene_id` | Ensembl gene identifier |
| `orf_id` | ORF identifier |
| `mut_id` | Mutation identifier |
| `ccsb_mutation_id` | CCSB mutation ID |
| `ccsb_allele_id` | CCSB allele ID |
| `ensembl_protein_id` | Ensembl protein identifier |

### Lab/Plate Information (Columns 17-48)
Information about clone locations and sequencing confirmation across different assay platforms:
- `collection` - Collection source
- `entry_plate_*`, `entry_well_*` - Entry clone plate/well locations
- `db_plate`, `db_well` - Destination plate/well
- `n2h_plate`, `n2h_well` - N2H assay plate/well
- `dualip_plate`, `dualip_well` - Dual IPA assay plate/well
- `mislocalization_plate`, `mislocalization_well` - Imaging assay plate/well
- `*_sequenced`, `*_sequence_confirmation_class` - Sequencing status and confirmation
- `*_seq_final` - Final sequencing status for each assay

### ClinVar Annotations (Columns 49-93)
Clinical significance annotations from ClinVar:
| Column | Description |
|--------|-------------|
| `#AlleleID` | ClinVar allele ID |
| `ClinicalSignificance` | Full clinical significance string |
| `ClinSigSimple` | Simplified clinical significance |
| `ReviewStatus` | ClinVar review status |
| `StarStatus` | ClinVar star rating |
| `clinvar_clnsig_clean` | Cleaned clinical significance categories |
| `clinvar_clnsig_clean_pp_strict` | Strict pathogenicity classification |
| `PhenotypeList` | Associated phenotypes |
| `Origin`, `OriginSimple` | Variant origin |

### Genomic Coordinates & Transcript Annotations (Columns 94-129)
| Column | Description |
|--------|-------------|
| `#chr`, `pos(1-based)` | GRCh38 coordinates |
| `hg19_chr`, `hg19_pos(1-based)` | GRCh37/hg19 coordinates |
| `aapos` | Amino acid position |
| `Ensembl_transcriptid` | Ensembl transcript ID |
| `Uniprot_acc`, `Uniprot_entry` | UniProt identifiers |
| `HGVSc_*`, `HGVSp_*` | HGVS nomenclature (snpEff and VEP) |
| `APPRIS`, `MANE`, `VEP_canonical` | Transcript annotation flags |
| `cds_strand` | Coding strand |
| `refcodon`, `codonpos` | Codon information |
| `Ancestral_allele` | Ancestral allele |
| `*Neandertal`, `Denisova` | Archaic human alleles |

### In Silico Pathogenicity Predictors (Columns 130-258)
Computational predictions of variant pathogenicity:

**Sequence-based predictors:**
- `SIFT_score`, `SIFT_pred` - SIFT predictions
- `SIFT4G_score`, `SIFT4G_pred` - SIFT4G predictions
- `Polyphen2_HDIV_*`, `Polyphen2_HVAR_*` - PolyPhen-2 predictions
- `MutationTaster_*` - MutationTaster predictions
- `MutationAssessor_*` - MutationAssessor predictions
- `PROVEAN_*` - PROVEAN predictions
- `VEST4_*` - VEST4 scores

**Ensemble/Meta predictors:**
- `MetaSVM_*`, `MetaLR_*`, `MetaRNN_*` - Meta-predictor scores
- `REVEL_score`, `REVEL_rankscore` - REVEL ensemble scores
- `M-CAP_*` - M-CAP predictions
- `MutPred_*` - MutPred predictions
- `MVP_*`, `gMVP_*` - MVP scores
- `MPC_*` - MPC scores
- `ClinPred_*` - ClinPred predictions
- `BayesDel_*` - BayesDel predictions

**Deep learning predictors:**
- `PrimateAI_*` - PrimateAI predictions
- `DEOGEN2_*` - DEOGEN2 predictions
- `ESM1b_*` - ESM-1b language model scores
- `AlphaMissense_*` - AlphaMissense predictions
- `PHACTboost_*`, `MutFormer_*`, `MutScore_*` - Additional ML predictors

**Conservation scores:**
- `CADD_raw`, `CADD_phred` - CADD scores
- `DANN_*` - DANN scores
- `fathmm-XF_*` - fathmm-XF predictions
- `Eigen-*` - Eigen scores
- `GERP++_*`, `GERP_91_mammals` - GERP conservation
- `phyloP*` - phyloP conservation scores (100way, 470way, 17way)
- `phastCons*` - phastCons conservation scores
- `bStatistic` - Background selection statistic

### Population Allele Frequencies (Columns 259-463)
Allele frequency data from multiple population databases:

**1000 Genomes Project (Phase 3):**
- `1000Gp3_AF` - Global allele frequency
- Population-specific: `_AFR_`, `_EUR_`, `_AMR_`, `_EAS_`, `_SAS_`

**TOPMed:**
- `TOPMed_frz8_AC`, `TOPMed_frz8_AN`, `TOPMed_frz8_AF`

**gnomAD v2.1.1 (exomes):**
- Control, non-neuro, and non-cancer subsets
- Population-specific frequencies (AFR, AMR, ASJ, EAS, FIN, NFE, SAS)
- `_AC` (allele count), `_AN` (allele number), `_AF` (frequency), `_nhomalt` (homozygotes)

**gnomAD v4.1 (joint):**
- Combined exome + genome data
- Additional populations: AMI (Amish), MID (Middle Eastern)

**ALFA (Allele Frequency Aggregator):**
- Population-specific frequencies from dbGaP studies

### ClinVar Metadata (Columns 464-471)
| Column | Description |
|--------|-------------|
| `clinvar_trait` | Associated traits/diseases |
| `clinvar_review` | Review status |
| `clinvar_hgvs` | HGVS nomenclature |
| `clinvar_MedGen_id` | MedGen identifier |
| `clinvar_OMIM_id` | OMIM identifier |
| `clinvar_Orphanet_id` | Orphanet identifier |
| `Interpro_domain` | InterPro domain annotations |

### Protein Structure Features (Columns 472-475)
| Column | Description |
|--------|-------------|
| `aaref_3`, `aaalt_3` | Reference and alternate amino acids (3-letter code) |
| `plddt` | AlphaFold pLDDT confidence score |
| `rsa` | Relative solvent accessibility |

### G2P DMS Assay Data (Columns 476-935)
Deep mutational scanning (DMS) data from various G2P assays. Each assay has two columns:
- `g2p_XXXXX-X-X avg` - Average functional score
- `g2p_XXXXX-X-X outliers⁺⁺` - Outlier classification

Assay IDs follow the pattern `g2p_NNNNN-letter-number` representing different experimental conditions and replicates.

### G2P Structural Annotations (Columns 936-1002)
Protein structure and post-translational modification annotations:

**Structure metrics:**
- `g2p_Accessible surface area (Å²)*`
- `g2p_AlphaFold confidence (pLDDT)`
- `g2p_Phi angle (degrees)*`, `g2p_Psi angle (degrees)*`
- `g2p_Secondary structure (DSSP 3-state)*`, `g2p_Secondary structure (DSSP 9-state)*`
- `g2p_Druggability score (fpocket)*`

**UniProt annotations:**
- Domain, region, motif, repeat information
- Active site, binding site, DNA binding annotations
- Transmembrane, signal peptide, transit peptide
- Disulfide bonds, cross-links

**Post-translational modifications (PTMs):**
- Phosphorylation, acetylation, methylation
- Glycosylation (O-GalNAc, O-GlcNAc)
- Ubiquitination, SUMOylation, lipidation
- Disease-associated PTMs, SNP-associated PTMs

**Molecular interactions:**
- Inter-chain and intra-chain interactions (PDB, AlphaFold2)
- Hydrogen bonds, salt bridges, disulfide bonds, non-bonded interactions

### AlphaFold Structural Metrics (Columns 1003-1006)
| Column | Description |
|--------|-------------|
| `plddt_f32` | pLDDT score (float32) |
| `rsa_f32` | RSA score (float32) |
| `pLDDT_Category` | Categorical pLDDT classification |
| `RSA_Category` | Categorical RSA classification (buried/exposed) |

### Disease & Inheritance Information (Columns 1007-1013)
| Column | Description |
|--------|-------------|
| `OMIM_IDs` | OMIM disease identifiers |
| `inheritance_pattern` | Mode of inheritance (AD, AR, XL, etc.) |
| `inheritance_source` | Source of inheritance annotation |
| `disease_modules` | Disease module classifications |
| `disease_module_count` | Number of associated disease modules |
| `OMIM_disease_names` | OMIM disease names |
| `splicing_variant` | Splicing variant flag |
