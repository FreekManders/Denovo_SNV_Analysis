#Requires vcftools
import re
import os
import pandas as pd
import datetime
import subprocess
from timeit import default_timer as timer
import argparse

###Start timer
start = timer()
time = datetime.datetime.now().strftime("%Y-%m-%d-%H-%M")
print "Start running script at: {0}".format(time)

###Parse arguments
parser = argparse.ArgumentParser(description = "This script identifies the snvs that are likely false positive and might therefore be manually inspected in IGV.")
parser.add_argument("-f", "--FILELIST", default = "/hpc/cog_bioinf/cuppen/project_data/Freek_SNV/FilelistFreek.txt", help = "The file list")
parser.add_argument("-o", "--OUTPUT_PATH", default = "/hpc/cog_bioinf/cuppen/project_data/Freek_SNV", help = "The output path")
parser.add_argument("--SNPSIFT", default = "/hpc/local/CentOS7/cog_bioinf/snpEff_v4.1h/SnpSift.jar", help = "The location of the snpsift jar.")
parser.add_argument("-g", "--GENOME", default = "/hpc/cog_bioinf/GENOMES/Homo_sapiens.GRCh37.GATK.illumina/Homo_sapiens.GRCh37.GATK.illumina.fasta", help = "The used genome in fasta format.")
parser.add_argument("--GATK", default = "/hpc/cog_bioinf/common_scripts/GenomeAnalysisTK-3.4-46/GenomeAnalysisTK.jar", help = "The GATK jar.")
args = parser.parse_args()

filelist = args.FILELIST
snpsift = args.SNPSIFT
gatk = args.GATK
genome = args.GENOME
folder = "{0}/Callableregion_NrMendelViols/TrueDenovo/SNVs".format(args.OUTPUT_PATH)
manualcheck = "{0}/Manualcheck".format(folder)

###Parse filelist
dffiles = pd.read_csv(filelist, sep = "\t", header = 0)
childrowsf = dffiles["Name"].str.contains("Child|child")
dffiles = dffiles[childrowsf]
dffiles["vcf"] = folder + "/" + dffiles["Sample"] + "_Truedenovo.vcf"
dffiles["table"] = folder + "/Manualcheck/" + dffiles["Sample"] + "temp.txt"
dffiles["vcfclose"] = folder + "/Manualcheck/" + dffiles["Sample"] + "_tempclosesites"
dffiles["callablevcfclose"] = dffiles["vcfclose"] + ".recode.vcf"
dffiles["tempcheckvcf"] = folder + "/Manualcheck/" + dffiles["Sample"] + "_tempmanualcheck.vcf"
dffiles["temp2checkvcf"] = folder + "/Manualcheck/" + dffiles["Sample"] + "_temp2manualcheck.vcf"
dffiles["temp3checkvcf"] = folder + "/Manualcheck/" + dffiles["Sample"] + "_temp3manualcheck.vcf"
dffiles["callabletemp3checkvcf"] = dffiles["temp3checkvcf"] + ".recode.vcf"
dffiles["checkvcf"] = folder + "/Manualcheck/" + dffiles["Sample"] + "_manualcheck"
dffiles["callablecheckvcf"] = dffiles["checkvcf"] + ".recode.vcf"


###For each sample check for sites that are likely false positive.
for index, row in dffiles.iterrows():
	sample = row["Sample"]
	vcf = row["vcf"]
	table = row["table"]
	tempcheckvcf = row["tempcheckvcf"]
	temp2checkvcf = row["temp2checkvcf"]
	temp3checkvcf = row["temp3checkvcf"]
	callabletemp3checkvcf = row["callabletemp3checkvcf"]
	vcfclose = row["vcfclose"]
	callablevcfclose = row["callablevcfclose"]
	checkvcf = row["checkvcf"]
	callablecheckvcf = row["callablecheckvcf"]

	###Ensure files with fps and tps exist
	fps = "{0}/FP/{1}_FPs.txt".format(folder, sample)
	open(fps, "a").close()
	tps = "{0}/TP/{1}_TPs.txt".format(folder, sample)
	open(tps, "a").close()

	if os.path.isfile(row["vcf"]):
		os.system("grep '^#' {0} > {1}".format(vcf, tempcheckvcf))
		
		###Grep the truedenovo vcf file to get sites for manual inspection, based on COSMIC or dbsnp ids or ambiguous phasing
		os.system("grep -v '^#' {0} | grep -e '^.*\t.*\tCOSM' -e '^.*\t.*\trs' -e 'PARENT=Ambiguous' >> {1}".format(vcf, tempcheckvcf))	
	
		###Find sites that are close together and add them to the manual inspection vcf
		#Create table from the truedenovo vcf
		command = 'java -Xmx4G -jar {0} extractFields -s "\t" -e "." {1} CHROM POS > {2}'.format(snpsift, vcf, table)
		os.system(command)

	
		#Select sites that are within 200bp from eachother
		df = pd.read_csv(table, sep = "\t", header = 0)
		df["CHROM2"] = df["#CHROM"].shift(-1)
		df["POS2"] = df["POS"].shift(-1)
		df["Diff"] = df["POS2"] - df["POS"]
		df["Close"] = (df["#CHROM"] == df["CHROM2"]) & (df["Diff"] <= 200)
		df["Close2"] = df["Close"].shift(1)
		df = df.loc[df.Close | df.Close2]
		df = df[["#CHROM", "POS"]]
		df.to_csv(table, sep = "\t", index = False)

		#Check if the number of sites within 200bp is not 0
		command = "grep -v '^#' {0} | wc -l".format(table)
		output = subprocess.check_output(command, shell=True)
		if output == "0\n":
			#No sites need to be added to the manual check sites
			os.rename(tempcheckvcf, temp2checkvcf)
			
		else:
			#Get the nearby sites from the truedenovo vcf
			command = "vcftools --positions {0} --vcf {1} --out {2} --recode --recode-INFO-all".format(table, vcf, vcfclose)
			os.system(command)
	
			#Add these sites to the vcf containing sites that need to be checked manually
			command = "java -Xmx4g -jar  {0} -T CombineVariants -R {1} --variant:foo {2} --variant:bar {3} -o {4} -genotypeMergeOptions PRIORITIZE -priority foo,bar".format(gatk, genome, tempcheckvcf, callablevcfclose, temp2checkvcf)
			os.system(command)
		
		###Remove known fps from sites to check manually
		command = "vcftools --exclude-positions {0} --vcf {1} --out {2} --recode --recode-INFO-all".format(fps, temp2checkvcf, temp3checkvcf)
		os.system(command)
		
		###Remove known tps from sites to check manually		
		command = "vcftools --exclude-positions {0} --vcf {1} --out {2} --recode --recode-INFO-all".format(tps, callabletemp3checkvcf, checkvcf)
		os.system(command)
				
		###The file with sites to check manually will only be retained if not empty
		command = "grep -v '^#' {0} | wc -l".format(callablecheckvcf)
		output = subprocess.check_output(command, shell=True)
		if output == "0\n":
			os.remove(callablecheckvcf)

###Remove temporary files
for fname in os.listdir(manualcheck):
	if re.search("temp|.log", fname):
		os.remove(os.path.join(manualcheck, fname))


###End timer
timeend = datetime.datetime.now().strftime("%Y-%m-%d-%H-%M")
print "Finished running script at: {0}".format(timeend)
end = timer()
execution_time = int(end - start)
execution_time = datetime.timedelta(seconds = execution_time)
print "Elapsed time: {0}".format(execution_time)
