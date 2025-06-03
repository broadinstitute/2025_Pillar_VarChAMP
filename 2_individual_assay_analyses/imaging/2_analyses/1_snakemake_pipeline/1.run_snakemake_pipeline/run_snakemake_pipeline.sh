#!/bin/bash

cd /home/shenrunx/igvf/varchamp/2025_Pillar_VarChAMP/2_individual_assay_analyses/imaging/2_analyses/1_snakemake_pipeline/2025_varchamp_snakemake/1.image_preprocess_qc/scripts
python 1_calc_plate_bg.py --batch_list "2025_01_27_Batch_13,2025_01_28_Batch_14" --input_dir "/home/shenrunx/igvf/varchamp/2025_Pillar_VarChAMP/2_individual_assay_analyses/imaging/1_inputs/1_snakemake_pipeline/cpg_imgs" --output_dir "../outputs/sum_stats_parquet" --workers 256

# cp snakemake_files/Snakefile_batch13 /home/shenrunx/igvf/varchamp/2025_Pillar_VarChAMP/2_individual_assay_analyses/imaging/2_analyses/1_snakemake_pipeline/2025_varchamp_snakemake/2.snakemake_pipeline/
# cp snakemake_files/Snakefile_batch14 /home/shenrunx/igvf/varchamp/2025_Pillar_VarChAMP/2_individual_assay_analyses/imaging/2_analyses/1_snakemake_pipeline/2025_varchamp_snakemake/2.snakemake_pipeline/
# cd /home/shenrunx/igvf/varchamp/2025_Pillar_VarChAMP/2_individual_assay_analyses/imaging/2_analyses/1_snakemake_pipeline/2025_varchamp_snakemake/2.snakemake_pipeline/

# ## Run batch 13
# snakemake \
#     --snakefile Snakefile_batch13 \
#     --directory /home/shenrunx/igvf/varchamp/2025_Pillar_VarChAMP/2_individual_assay_analyses/imaging/2_analyses/1_snakemake_pipeline/2025_varchamp_snakemake/2.snakemake_pipeline/ \
#     --cores 256 &> /home/shenrunx/igvf/varchamp/2025_Pillar_VarChAMP/2_individual_assay_analyses/imaging/3_outputs/1_snakemake_pipeline/pipeline_outputs/snakemake_logs/snakemake_batch13.log

## Run batch 14
# snakemake \
#     --snakefile Snakefile_batch14 \
#     --cores 256 &> /home/shenrunx/igvf/varchamp/2025_Pillar_VarChAMP/2_individual_assay_analyses/imaging/3_outputs/1_snakemake_pipeline/pipeline_outputs/snakemake_logs/snakemake_batch14.log