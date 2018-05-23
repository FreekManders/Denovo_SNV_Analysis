import os
import pandas as pd
import datetime
import argparse
from timeit import default_timer as timer

###Start timer
start = timer()
time = datetime.datetime.now().strftime("%Y-%m-%d-%H-%M")
print "Start running script at: {0}".format(time)

###Parse arguments
parser = argparse.ArgumentParser(description = "This script filters vcf files on 'callable' loci and then runs it through PhaseByTransmission.")
parser.add_argument("-f", "--FILELIST", default = "/hpc/cog_bioinf/cuppen/project_data/Freek_SNV/FilelistFreek.txt", help = "The file list")
parser.add_argument("-o", "--OUTPUT_PATH", default = "/hpc/cog_bioinf/cuppen/project_data/Freek_SNV", help = "The output path")
args = parser.parse_args()

filelist = args.FILELIST
folder = "{0}/Callableregion_NrMendelViols/TrueDenovo/SNVs".format(args.OUTPUT_PATH)

###Parse filelist
dffiles = pd.read_csv(filelist, sep = "\t", header = 0)
childrowsf = dffiles["Name"].str.contains("Child|child")
dffiles = dffiles[childrowsf]

###Remove false positives
for sample in dffiles["Sample"]:
	vcf = "{0}/{1}_Truedenovo.vcf".format(folder, sample)
	fps = "{0}/FP/{1}_FPs.txt".format(folder, sample)
	tempvcf = "{0}/FPfilterlogs/{1}_filterFPs_{2}".format(folder, sample, time)
	callabletempvcf = "{0}.recode.vcf".format(tempvcf)
	if os.path.isfile(vcf) and os.path.isfile(fps):
		command = "vcftools --exclude-positions {0} --vcf {1} --out {2} --recode --recode-INFO-all".format(fps, vcf, tempvcf)
		os.system(command)
		os.remove(vcf)
		os.rename(callabletempvcf, vcf)
		print "Removed false positives for sample: {0}".format(sample)
		
###End timer
timeend = datetime.datetime.now().strftime("%Y-%m-%d-%H-%M")
print "Finished running script at: {0}".format(timeend)
end = timer()
execution_time = int(end - start)
execution_time = datetime.timedelta(seconds = execution_time)
print "Elapsed time: {0}".format(execution_time)

