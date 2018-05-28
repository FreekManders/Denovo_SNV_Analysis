import os
import pandas as pd
from multiprocessing import Pool as ThreadPool
import argparse
import datetime
from timeit import default_timer as timer

###Start timer
start = timer()
time = datetime.datetime.now().strftime("%Y-%m-%d-%H-%M")
print "Start running script at: {0}".format(time)

###Parse arguments
parser = argparse.ArgumentParser(description = "This script annotates a vcf with DANN scores.")
parser.add_argument("-f", "--FILELIST", default = "/hpc/cog_bioinf/cuppen/project_data/Freek_SNV/FilelistFreek.txt", help = "The file list")
parser.add_argument("-o", "--OUTPUT_PATH", default = "/hpc/cog_bioinf/cuppen/project_data/Freek_SNV", help = "The output path")
parser.add_argument("-w", "--OVERWRITE", default = "False", help = "Overwrite the results if they already exist.")
parser.add_argument("--DANN_FOLDER", default = "/hpc/cog_bioinf/cuppen/project_data/Freek_SNV/Downloaded_software/DANN", help = "The location of the DANN files.")
parser.add_argument("-t", "--THREADS", default = 1, type = int, help = "Number of threads to use")
args = parser.parse_args()

filelist = args.FILELIST
folder = "{0}/Callableregion_NrMendelViols/TrueDenovo/SNVs".format(args.OUTPUT_PATH)
DANN = "{0}/DANN_whole_genome_SNVs.tsv.bgz".format(args.DANN_FOLDER)
DANNhdr = "{0}/DANNannotation.hdr".format(args.DANN_FOLDER)
annotatedf = "{0}/Tables/DANNannotated.txt".format(folder)

###Parse filelist
dffiles = pd.read_csv(filelist, sep = "\t", header = 0)
childrowsf = dffiles["Name"].str.contains("Child|child")
dffiles = dffiles[childrowsf]
dffiles["vcf"] = folder + "/" + dffiles["Sample"] + "_Truedenovo.vcf"
dffiles["vcftempout"] = folder + "/Temp/" + dffiles["Sample"] + "_Truedenovo.vcf"

###Check if annotation files exist
if not os.path.isfile(DANN) or not os.path.isfile(DANNhdr):
	raise IOError("DANN file or the annotation.hdr file does not exist")

###Remove samples that have already been run
if not args.OVERWRITE == "true" and os.path.isfile(annotatedf):
	annotated = pd.read_csv(annotatedf, header = None)
	alreadyrun = dffiles["Sample"].isin(annotated[0])
	dffiles = dffiles[~alreadyrun]

###Put filepaths in a nested list. This list conains a list with filepaths for each sample.
lfiles = dffiles[["Sample", "vcf", "vcftempout"]].values.tolist()

###Annotation function
def AnnotateDANN(lfiles_sample):
	sample, vcf, vcftempout = lfiles_sample
	if os.path.isfile(vcf):
		command = "bcftools annotate -a {0} -c CHROM,POS,REF,ALT,DANN {1} -h {2} -o {3}".format(DANN, vcf, DANNhdr, vcftempout)
		os.system(command)
		if os.stat(vcftempout).st_size == 0:
			os.remove(vcftempout)
			raise Exception("Temporary annotated vcf is empty. Something likely went wrong with bcftools annotate.")
		os.remove(vcf)
		os.rename(vcftempout, vcf)
		print "Finished annotating sample: {0}".format(sample)
		with open(annotatedf, "a") as annotated:
			annotated.write("{0}\n".format(sample))


#Map function over the nested list
pool = ThreadPool(args.THREADS)
results = pool.map(AnnotateDANN, lfiles)
pool.close()
pool.join()


#Removes duplicates from the list of already annotated samples.
with open(annotatedf, "r") as annotated:
	uniqlines = sorted(set(annotated.readlines()))

with open(annotatedf, "w") as annotated:
	annotated.writelines(uniqlines)

###End timer
timeend = datetime.datetime.now().strftime("%Y-%m-%d-%H-%M")
print "Finished running script at: {0}".format(timeend)
end = timer()
execution_time = int(end - start)
execution_time = datetime.timedelta(seconds = execution_time)
print "Elapsed time: {0}".format(execution_time)
