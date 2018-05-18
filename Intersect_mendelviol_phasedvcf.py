import os
import argparse
import datetime
from timeit import default_timer as timer

###Start timer
start = timer()
time = datetime.datetime.now().strftime("%Y-%m-%d-%H-%M")
print "Start running script at: {0}".format(time)

###Parse arguments
parser = argparse.ArgumentParser(description = "This script intersects a mendelviol file with a phased vcf.")
parser.add_argument("-f", "--FILELIST", default = "/hpc/cog_bioinf/cuppen/project_data/Freek_SNV/FilelistFreek.txt", help = "The file list")
parser.add_argument("-o", "--OUTPUT_PATH", default = "/hpc/cog_bioinf/cuppen/project_data/Freek_SNV", help = "The output path")
parser.add_argument("-w", "--OVERWRITE", default = "False", help = "Overwrite the results if they already exist.")
args = parser.parse_args()

filelist = args.FILELIST
outdir = "{0}/Callableregion_NrMendelViols".format(args.OUTPUT_PATH)
		
###Remove the missing GTs from the MendelViol files.
with open(filelist) as f1:
	next(f1)
	for line in f1:
		line = line.rstrip()
		columns = line.split("\t")
		sample = columns[0]
		family = columns[3]
		relation = columns[6]
		relation = relation.lower()
		if "child" in relation:
			mendel_fname = "{0}/PBT/{1}_fullgenotypes.MendelViol".format(outdir, family)
			vcf_fname = "{0}/PBT/{1}.phased.vcf.gz".format(outdir, family)
			denovovcf = "{0}/Denovo/{1}_denovo".format(outdir, sample)
			denovovcf_fname = "{0}.recode.vcf".format(denovovcf)
			if os.path.isfile(mendel_fname) and os.path.isfile(vcf_fname) and (not os.path.isfile(denovovcf_fname) or args.OVERWRITE == "true"):
				command = "vcftools --gzvcf {0} --positions {1} --out {2} --recode --recode-INFO-all".format(vcf_fname, mendel_fname, denovovcf)
				os.system(command)
	
###End timer
timeend = datetime.datetime.now().strftime("%Y-%m-%d-%H-%M")
print "Finished running script at: {0}".format(timeend)
end = timer()
execution_time = int(end - start)
execution_time = datetime.timedelta(seconds = execution_time)
print "Elapsed time: {0}".format(execution_time)
