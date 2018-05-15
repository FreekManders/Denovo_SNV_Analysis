#!/bin/bash

#$ -S /bin/bash
#$ -l h_rt=00:20:00
#$ -l h_vmem=10G
#$ -cwd
#$ -o /hpc/cog_bioinf/cuppen/project_data/Freek_SNV/pipeline_output
#$ -e /hpc/cog_bioinf/cuppen/project_data/Freek_SNV/pipeline_errors
#$ -M fmanders@umcutrecht.nl
#$ -m beas

echo "Running FilterDenovo on: $(date)"

module load python/2.7.10

python /hpc/cog_bioinf/cuppen/project_data/Freek_SNV/pipeline/Denovo_snv_analysis.py

module unload python/2.7.10

