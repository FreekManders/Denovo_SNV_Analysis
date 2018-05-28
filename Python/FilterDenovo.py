import os
import pandas as pd
import argparse
import datetime
from timeit import default_timer as timer

"""#Filters:
Check whether everything passed the previous filters.
Check whether genotype quality of all samples is above 50.
Check whether the read depth of all samples is above 10 and below 120.
Check whether the child is heterogenous.
Check whether the Allele depth of both child alleles is at least 5.
Check whether the Allele frequency of the alternate allele is at least 0.3.
Check whether the Allele frequency of the reference allele is at least 0.2.
Check whether the parents have no informative reads on the alternative allele.
"""

###Start timer
start = timer()
time = datetime.datetime.now().strftime("%Y-%m-%d-%H-%M")
print "Start running script at: {0}".format(time)

###Parse arguments
parser = argparse.ArgumentParser(description = "This script extracts the loci that were 'Callable' in all members in a parent-offspring trio These loci are extracted from bed files generated by GATKs CallableLoci.")
parser.add_argument("-f", "--FILELIST", default = "/hpc/cog_bioinf/cuppen/project_data/Freek_SNV/FilelistFreek.txt", help = "The file list")
parser.add_argument("-o", "--OUTPUT_PATH", default = "/hpc/cog_bioinf/cuppen/project_data/Freek_SNV", help = "The output path")
parser.add_argument("-w", "--OVERWRITE", default = "False", help = "Overwrite the results if they already exist.")
parser.add_argument("--SNPSIFT", default = "/hpc/local/CentOS7/cog_bioinf/snpEff_v4.1h/SnpSift.jar", help = "The location of the snpsift jar.")
args = parser.parse_args()

filelist = args.FILELIST
outdir = "{0}/Callableregion_NrMendelViols".format(args.OUTPUT_PATH)
snpsift = args.SNPSIFT

###The filter settings
filtersautosome = "(FILTER = 'PASS') & (GEN[?].GQ > 50) & (GEN[?].DP > 10) & (GEN[?].DP < 120) & isHet( GEN[2] ) & (GEN[2].AD[0] >= 5) & (GEN[2].AD[1] >= 5) & ((GEN[2].AD[0]/(GEN[2].AD[0] + GEN[2].AD[1] + 0.01)) > 0.2) & ((GEN[2].AD[1]/(GEN[2].AD[0] + GEN[2].AD[1] + 0.01)) > 0.3) & (GEN[0].AD[1] = 0) & (GEN[1].AD[1] = 0) & !(CHROM = 'X') & !(CHROM = 'Y') & !(CHROM = 'MT')"

###Get readdepth from wgs quality control tables
wgs = pd.DataFrame()
runs = []
with open(filelist) as f1:
	f1.next()
	for line in f1:
		line = line.rstrip()
		columns = line.split("\t")
		run = columns[1]	
		if run not in runs:
			runs.append(run)
			wgsf = "{0}/QCStats/WGSMetrics_summary.txt".format(run)
			wgsrun = pd.read_csv(wgsf, header = 0, sep = "\t")
			wgs = wgs.append(wgsrun)
wgs["sample "].replace("-DNA-1_dedup", "", inplace = True, regex = True)
wgs["family"] = wgs["sample "].replace("-.*", "", inplace = False, regex = True)
wgs["relation"] = wgs["sample "].replace(".*-", "", inplace = False, regex = True).apply(lambda x: x[0])
wgs["maxdp"] = wgs["MEDIAN_COVERAGE"] + 2 * wgs["SD_COVERAGE"]

maxdps = []
###Filter for each child sample
childnr = 1
with open(filelist) as f1:
	next(f1)
	for line in f1:
		line = line.rstrip()
		columns = line.split("\t")
		sample = columns[0]
		family = columns[3]
		sex = columns[4]
		relation = columns[6]
		relation = relation.lower()
		if "child" in relation:
			denovovcf = "{0}/Denovo/{1}_denovo.recode.vcf".format(outdir, sample)
			print denovovcf
			
			#Keep only SNVs
			SNVdenovovcf = "{0}/Denovo/SNVs/{1}_denovo_snv".format(outdir, sample)
			callable_SNVdenovovcf = "{0}.recode.vcf".format(SNVdenovovcf)
			if os.path.isfile(denovovcf) and (not os.path.isfile(callable_SNVdenovovcf) or args.OVERWRITE == "true"):
				command = "vcftools --vcf {0} --remove-indels --out {1} --recode --recode-INFO-all".format(denovovcf, SNVdenovovcf)
				os.system(command)
				print "Created a vcf with only SNVs for sample: {0}".format(sample)
				
			#Determine maximum dps for the filter
			childwgs = wgs.loc[wgs["sample "] == sample]
			childmaxdp = childwgs.iloc[0]["maxdp"]
			fatherwgs = wgs.loc[(wgs["family"] == family) & (wgs["relation"] == "1")]
			fathermaxdp = fatherwgs.iloc[0]["maxdp"]
			motherwgs = wgs.loc[(wgs["family"] == family) & (wgs["relation"] == "2")]
			mothermaxdp = motherwgs.iloc[0]["maxdp"]
			filtersautosome2 = "{0} & (GEN[0].DP <= {1}) & (GEN[1].DP <= {2}) & (GEN[2].DP <= {3})".format(filtersautosome, fathermaxdp, mothermaxdp, childmaxdp)
			maxdps.append([sample, fathermaxdp, mothermaxdp, childmaxdp])

			#Filter the denovo snv vcf file.
			truedenovovcfdpf = "{0}/TrueDenovo/SNVs/{1}_Truedenovo.vcf".format(outdir, sample)
			if os.path.isfile(callable_SNVdenovovcf) and (not os.path.isfile(truedenovovcfdpf) or args.OVERWRITE == "true"):
				command = 'java -Xmx4G -jar {0} filter "{1}" {2} > {3}'.format(snpsift, filtersautosome2, callable_SNVdenovovcf, truedenovovcfdpf)
				os.system(command)
				print "Filtered the denovo file for sample {0}".format(sample)
			
			#Keep only indels
			Indeldenovovcf = "{0}/Denovo/Indels/{1}_denovo_indel".format(outdir, sample)
			callable_Indeldenovovcf = "{0}.recode.vcf".format(Indeldenovovcf)
			if os.path.isfile(denovovcf) and (not os.path.isfile(callable_Indeldenovovcf) or args.OVERWRITE == "true"):
				command = "vcftools --vcf {0} --keep-only-indels --out {1} --recode --recode-INFO-all".format(denovovcf, Indeldenovovcf)
				os.system(command)
				print "Created a vcf with only the indels for sample {0}".format(sample)

###Create a table of the maximum DPs
maxdpf = "{0}/TrueDenovo/SNVs/Tables/maxalloweddepths.txt".format(outdir)
maxdps = pd.DataFrame(maxdps, columns = ["sample", "fathermaxdp", "mothermaxdp", "childmaxdp"])
maxdps.to_csv(maxdpf, header = True, index = False, sep = "\t")
print "Created a table with the maximum DPs"
	
###End timer
timeend = datetime.datetime.now().strftime("%Y-%m-%d-%H-%M")
print "Finished running script at: {0}".format(timeend)
end = timer()
execution_time = int(end - start)
execution_time = datetime.timedelta(seconds = execution_time)
print "Elapsed time: {0}".format(execution_time)