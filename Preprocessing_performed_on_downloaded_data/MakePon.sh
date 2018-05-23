#!/bin/bash

#$ -S /bin/bash
#$ -l h_rt=04:30:00
#$ -l h_vmem=10G
#$ -pe threaded 1
#$ -cwd
#$ -o /hpc/cog_bioinf/cuppen/project_data/Freek_SNV/Callableregion_NrMendelViols/Logs/MakePonfiles_output
#$ -e /hpc/cog_bioinf/cuppen/project_data/Freek_SNV/Callableregion_NrMendelViols/Logs/MakePonfiles_errors
#$ -M fmanders@umcutrecht.nl
#$ -m beas

echo "Running MakePon.sh on: $(date)"

module load tabix/0.2.6
module load vcfbcf/vcftools/c7a7337



#Remove indels and sites with filter not PASS
vcftools --gzvcf /hpc/cog_bioinf/common_dbs/HMF_PON/PON_v2.0.vcf.gz --remove-indels --remove-filtered-all --stdout --recode --recode-INFO-all | bgzip > /hpc/cog_bioinf/cuppen/project_data/Freek_SNV/Downloaded_software/HMF-PON/PON_v2.0_pass_snv.vcf.gz
tabix -p vcf /hpc/cog_bioinf/cuppen/project_data/Freek_SNV/Downloaded_software/HMF-PON/PON_v2.0_pass_snv.vcf.gz


module unload tabix/0.2.6
module unload vcfbcf/vcftools/c7a7337

