#!/bin/bash

cp snakemake_files/Snakefile_batch13 /home/shenrunx/igvf/varchamp/2025_Pillar_VarChAMP/2_individual_assay_analyses/imaging/2_analyses/1_snakemake_pipeline/2025_varchamp_snakemake/2.snakemake_pipeline/
cp snakemake_files/Snakefile_batch14 /home/shenrunx/igvf/varchamp/2025_Pillar_VarChAMP/2_individual_assay_analyses/imaging/2_analyses/1_snakemake_pipeline/2025_varchamp_snakemake/2.snakemake_pipeline/
cd /home/shenrunx/igvf/varchamp/2025_Pillar_VarChAMP/2_individual_assay_analyses/imaging/2_analyses/1_snakemake_pipeline/2025_varchamp_snakemake/2.snakemake_pipeline/

## Run batch 13
snakemake \
    --snakefile Snakefile_batch13 \
    --directory /home/shenrunx/igvf/varchamp/2025_Pillar_VarChAMP/2_individual_assay_analyses/imaging/2_analyses/1_snakemake_pipeline/2025_varchamp_snakemake/2.snakemake_pipeline/ \
    --cores 256 &> /home/shenrunx/igvf/varchamp/2025_Pillar_VarChAMP/2_individual_assay_analyses/imaging/3_outputs/1_snakemake_pipeline/pipeline_outputs/snakemake_logs/snakemake_batch13.log

## Run batch 14
snakemake \
    --snakefile Snakefile_batch14 \
    --cores 256 &> /home/shenrunx/igvf/varchamp/2025_Pillar_VarChAMP/2_individual_assay_analyses/imaging/3_outputs/1_snakemake_pipeline/pipeline_outputs/snakemake_logs/snakemake_batch14.log