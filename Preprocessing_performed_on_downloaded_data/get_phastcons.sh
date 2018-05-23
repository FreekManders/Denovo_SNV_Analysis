#!/bin/bash

#$ -S /bin/bash
#$ -l h_rt=00:05:00
#$ -l h_vmem=10G
#$ -pe threaded 1
#$ -cwd
#$ -o /hpc/cog_bioinf/cuppen/project_data/Freek_SNV/Downloaded_software/Logs/phastCons46way_output
#$ -e /hpc/cog_bioinf/cuppen/project_data/Freek_SNV/Downloaded_software/Logs/phastCons46way_errors
#$ -M fmanders@umcutrecht.nl
#$ -m beas


cd /hpc/cog_bioinf/cuppen/project_data/Freek_SNV/



#The following repositories and versions were used:
# * GNU Guix upstream:         d000264299864480aa8b80357a3c6ae5ff2a5e74
# * UMCU additional packages:  1c49c8b2d3b9c33def43a30f53c73941317ffd7b
#guixr package -i coreutils grep bedops htslib -p guix-profile

#The source should work if submitting this script as a job. Otherwise use guixr load profile and then run the script as a bash file. In this case: Do Not submit as a job.
#guixr load-profile guix-profile --<<EOF
#source /hpc/cog_bioinf/cuppen/project_data/Freek_SNV/guix-profile/etc/profile

cd Downloaded_software/phastCons46way/

conservation="phastCons46way"

for i in {1..22}
do
file="chr${i}.${conservation}"
rsync -avz --progress rsync://hgdownload.cse.ucsc.edu/goldenPath/hg19/${conservation}/vertebrate/${file}.wigFix.gz ${file}.wigFix.gz
gunzip ${file}.wigFix.gz
echo "Unzipped the file for chromosome ${i}"
wig2bed --do-not-sort < ${file}.wigFix > ${file}.bed
echo "Converted the wiggle file for chromosome ${i} to a bed file"
rm ${file}.wigFix
awk -F "\t" 'BEGIN {OFS = FS} {gsub(/^chr/,""); print $1, $2, $3, $5}' ${file}.bed >> ${conservation}.bed #remove the "chr" in the chromosome column, remove column 4 and append all chromosomes.
echo "Finished with the file for chromosome ${i}"
done

bgzip ${conservation}.bed
echo "Finished compressing the ${conservation} file"
tabix -p bed ${conservation}.bed.gz
echo "Finished indexing the ${conservation} file"
#EOF
