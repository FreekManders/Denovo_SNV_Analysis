import os
import pandas as pd
import datetime
import argparse
from timeit import default_timer as timer
import subprocess
import signal

###Start timer
start = timer()
time = datetime.datetime.now().strftime("%Y-%m-%d-%H-%M")
print "Start running script at: {0}".format(time)

###Parse arguments
parser = argparse.ArgumentParser(description = "For each sample this script reads the vcf file with true denovo mutations. It then writes two vcf files containing the mutations linked to the father and mother respectively. These vcfs can be used in mutationalpatterns.")
parser.add_argument("-f", "--FILELIST", default = "/hpc/cog_bioinf/cuppen/project_data/Freek_SNV/FilelistFreek.txt", help = "The file list")
parser.add_argument("-o", "--OUTPUT_PATH", default = "/hpc/cog_bioinf/cuppen/project_data/Freek_SNV", help = "The output path")
parser.add_argument("-w", "--OVERWRITE", default = "False", help = "Overwrite the results if they already exist.")
parser.add_argument("--SNPSIFT", default = "/hpc/local/CentOS7/cog_bioinf/snpEff_v4.1h/SnpSift.jar", help = "The location of the snpsift jar.")
args = parser.parse_args()

filelist = args.FILELIST
folder = "{0}/Callableregion_NrMendelViols/TrueDenovo/SNVs".format(args.OUTPUT_PATH)
snpsift = args.SNPSIFT

###Read the filelist
dffiles = pd.read_csv(filelist, sep = "\t", header = 0)
childrowsf = dffiles["Name"].str.contains("Child|child")
dffiles = dffiles[childrowsf]
dffiles["vcf"] = folder + "/" + dffiles["Sample"] + "_Truedenovo.vcf"
dffiles["vcffather"] = folder + "/Parentaloriginmuts/" + dffiles["Sample"] + "_Truedenovofather.vcf"
dffiles["vcfmother"] = folder + "/Parentaloriginmuts/" + dffiles["Sample"] + "_Truedenovomother.vcf"

###Filter snvs with a known parent of origin.
for index, row in dffiles.iterrows():
	sample, vcf, bam, vcfmother, vcffather = row[["Sample", "vcf", "bam", "vcfmother", "vcffather"]].values.tolist()
	if os.path.isfile(vcf):
		with open(vcf) as vcf_obj:#check if vcf is phased
			phased = "PARENT" in vcf_obj.read()
		
		if phased and (not os.path.isfile(vcffather) or args.OVERWRITE == "true"):
			command = """java -Xmx4G -jar {0} filter "( PARENT has 'Father' )" {1} > {2}""".format(snpsift, vcf, vcffather)
			os.system(command)
			print "Filtered paternal muts for sample: {0}".format(sample)
		if phased and (not os.path.isfile(vcfmother) or args.OVERWRITE == "true"):		
			command = """java -Xmx4G -jar {0} filter "( PARENT has 'Mother' )" {1} > {2}""".format(snpsift, vcf, vcfmother)
			os.system(command)
			print "Filtered maternal muts for sample: {0}".format(sample)
		

###End timer
timeend = datetime.datetime.now().strftime("%Y-%m-%d-%H-%M")
print "Finished running script at: {0}".format(timeend)
end = timer()
execution_time = int(end - start)
execution_time = datetime.timedelta(seconds = execution_time)
print "Elapsed time: {0}".format(execution_time)


