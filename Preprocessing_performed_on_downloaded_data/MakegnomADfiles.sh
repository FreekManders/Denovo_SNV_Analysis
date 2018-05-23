#!/bin/bash

#$ -S /bin/bash
#$ -l h_rt=00:30:00
#$ -l h_vmem=10G
#$ -pe threaded 1
#$ -cwd
#$ -o /hpc/cog_bioinf/cuppen/project_data/Freek_SNV/Callableregion_NrMendelViols/Logs/MakegnomADfiles_output
#$ -e /hpc/cog_bioinf/cuppen/project_data/Freek_SNV/Callableregion_NrMendelViols/Logs/MakegnomADfiles_errors
#$ -M fmanders@umcutrecht.nl
#$ -m beas

echo "Running MakegnomAD.sh on: $(date)"

module load tabix/0.2.6
module load vcfbcf/vcftools/c7a7337
module load python

#Create vcf with only filter == pass sites
#java -Xmx4G -jar /hpc/local/CentOS7/cog_bioinf/snpEff_v4.1h/SnpSift.jar filter "(FILTER = 'PASS')" /hpc/cog_bioinf/common_dbs/GNOMAD/v2.0.2/gnomad.genomes.r2.0.2.sites.vcf.gz > /hpc/cog_bioinf/cuppen/project_data/Freek_SNV/Downloaded_software/gnomAD/gnomADpass.vcf
#bgzip /hpc/cog_bioinf/cuppen/project_data/Freek_SNV/Downloaded_software/gnomAD/gnomADpass.vcf
#tabix -p vcf /hpc/cog_bioinf/cuppen/project_data/Freek_SNV/Downloaded_software/gnomAD/gnomADpass.vcf.gz

#Remove indels
#vcftools --gzvcf /hpc/cog_bioinf/cuppen/project_data/Freek_SNV/Downloaded_software/gnomAD/gnomADpass.vcf.gz --remove-indels --stdout --recode --recode-INFO-all | bgzip > /hpc/cog_bioinf/cuppen/project_data/Freek_SNV/Downloaded_software/gnomAD/gnomADpass_snv.recode.vcf.gz
#tabix -p vcf /hpc/cog_bioinf/cuppen/project_data/Freek_SNV/Downloaded_software/gnomAD/gnomADpass_snv.recode.vcf.gz

#python /hpc/cog_bioinf/cuppen/project_data/Freek_SNV/Downloaded_software/gnomAD/makegnomADfiles2.py

module unload tabix/0.2.6
module unload vcfbcf/vcftools/c7a7337
module unload python
