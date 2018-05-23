import os
import pandas as pd
import datetime
from timeit import default_timer as timer


#Start timer
start = timer()
time = datetime.datetime.now().strftime("%Y-%m-%d-%H-%M")
print "Start running script at: {0}".format(time)

#Create the file paths.
snpsift = "/hpc/local/CentOS7/cog_bioinf/snpEff_v4.1h/SnpSift.jar" 
gnomADfolder = "/hpc/cog_bioinf/cuppen/project_data/Freek_SNV/Downloaded_software/gnomAD"
gnomAD = "/hpc/cog_bioinf/cuppen/project_data/Freek_SNV/Downloaded_software/gnomAD/gnomADpass_snv.recode.vcf.gz"
gatk = "/hpc/cog_bioinf/common_scripts/GenomeAnalysisTK-3.4-46/GenomeAnalysisTK.jar"
genome = "/hpc/cog_bioinf/GENOMES/Homo_sapiens.GRCh37.GATK.illumina/Homo_sapiens.GRCh37.GATK.illumina.fasta"

bashsettings = "#!/bin/bash\n#$ -S /bin/bash\n#$ -l h_rt=00:45:00\n#$ -l h_vmem=7G\n#$ -cwd\n#$ -N MakegnomADfiles\n#$ -o /hpc/cog_bioinf/cuppen/project_data/Freek_SNV/Callableregion_NrMendelViols/Logs/MakegnomADfiles_output\n#$ -e /hpc/cog_bioinf/cuppen/project_data/Freek_SNV/Callableregion_NrMendelViols/Logs/MakegnomADfiles_errors\n#$ -M fmanders@umcutrecht.nl\n#$ -m beas\nmodule load tabix/0.2.6"


#Per chromosome
for i in range(1, 23):

	gnomadchrom = "{0}/gnomADpass_snv.recode.{1}.vcf".format(gnomADfolder, i)
	gnomadchromtable = "{0}/gnomADpass_snv_table_{1}.txt".format(gnomADfolder, i)
	
	#Split vcf
	command1 = "tabix -h {0} {1} > {2}".format(gnomAD, i, gnomadchrom)

	#Extract to table
	command2 = """java -Xmx4g -jar {0} -T VariantsToTable -R {1} -V {2} -F CHROM -F POS -F REF -F ALT -F AC -F AF -F AN -o {3} --splitMultiAllelic""".format(gatk, genome, gnomadchrom, gnomadchromtable)

	#Compress and remove intermediate files
	command3 = "bgzip {0}".format(gnomadchromtable)
	command4 = "tabix -s1 -b2 -e2 {0}.gz".format(gnomadchromtable)
	command5 = "rm {0}".format(gnomadchrom)
	command6 = "rm {0}.idx".format(gnomadchrom)
	
	#Create the bash script
	bashname = "{0}/Jobs/chr{1}_{2}.sh".format(gnomADfolder, i, time)
	fullbash = "{0}\n{1}\n{2}\n{3}\n{4}\n{5}\n{6}\nmodule unload tabix/0.2.6".format(bashsettings, command1, command2, command3, command4, command5, command6)
	with open(bashname, "w") as bash:
		bash.write(fullbash)
	os.system("qsub {0}".format(bashname))




	
#End timer		
timeend = datetime.datetime.now().strftime("%Y-%m-%d-%H-%M")
print "Finished running script at: {0}".format(timeend)
end = timer()
execution_time = int(end - start)
execution_time = datetime.timedelta(seconds = execution_time)
print "Elapsed time: {0}".format(execution_time)
