import os
import subprocess
import signal
import argparse
import datetime
from timeit import default_timer as timer


###Start timer
start = timer()
time = datetime.datetime.now().strftime("%Y-%m-%d-%H-%M")
print "Start running script at: {0}".format(time)

###Parse arguments
parser = argparse.ArgumentParser(description = "This script is meant for families with three children. For each child it removes the samples from the other kids from a vcf. It filters the vcf on callable loci. It removes the other kids from a ped file. It runs this filtered vcf through PhaseByTransmission. It then removes the missing genotypes and intersect the vcf with the mendelviol file.")
parser.add_argument("-f", "--FILELIST", default = "/hpc/cog_bioinf/cuppen/project_data/Freek_SNV/FilelistFreek.txt", help = "The file list")
parser.add_argument("-o", "--OUTPUT_PATH", default = "/hpc/cog_bioinf/cuppen/project_data/Freek_SNV", help = "The output path")
parser.add_argument("-w", "--OVERWRITE", default = "False", help = "Overwrite the results if they already exist.")
parser.add_argument("-g", "--GENOME", default = "/hpc/cog_bioinf/GENOMES/Homo_sapiens.GRCh37.GATK.illumina/Homo_sapiens.GRCh37.GATK.illumina.fasta", help = "The used genome in fasta format.")
parser.add_argument("--GATK", default = "/hpc/cog_bioinf/common_scripts/GenomeAnalysisTK-3.4-46/GenomeAnalysisTK.jar", help = "The GATK jar.")
parser.add_argument("--MUTATIONPRIOR", default = "1.0E-4", help = "The mutation prior that is being used.")
parser.add_argument("--FAMILY", default = "MP14", help = "Family to run this script on")
args = parser.parse_args()

filelist = args.FILELIST
genome = args.GENOME
jar = args.GATK
family = args.FAMILY
outdir = args.OUTPUT_PATH
famdir = "{0}/Callableregion_NrMendelViols/{1}".format(outdir, family)

###Find the correct vcf file
with open(filelist) as f1:
	next(f1)
	for line in f1:
		line = line.rstrip()
		columns = line.split("\t")
		if columns[3] == family:
			run = columns[1]
			if os.path.isdir("{0}/vcf_new_filter".format(run)):
				vcfdir = "vcf_new_filter"
			else:
				vcfdir = "vcf"
			vcffile = "{0}/{1}/{2}/{2}.filtered_variants_snpEff_snpSift_Cosmicv76_GoNLv5.vcf.gz".format(run, vcfdir, family)
			if not os.path.isfile(vcffile):
				vcffile = "{0}/{1}/{2}/{2}.filtered_variants_snpEff_snpSift_Cosmicv76_GoNLv5.vcf".format(run, vcfdir, family)
				if not os.path.isfile(vcffile):
					print "vcf file does not exists."
			break

###Create family directory
if not os.path.exists(famdir):
	os.mkdir(famdir)

###Find the sample ids in the vcf
sampleids = subprocess.check_output(["vcf-query", "-l", vcffile], preexec_fn=lambda:signal.signal(signal.SIGPIPE, signal.SIG_DFL))
sampleids = sampleids.split("\n")
del sampleids[-1] #Remove empty newline

###Main body of script
idstoremove = [[3, 4], [2, 4], [2, 3]]
for childnr in [1,2,3]:
	sample = sampleids[childnr+1].rstrip("-DNA-1")
	#Remove the other samples from the vcf
	idstoremove_child = idstoremove[childnr-1]
	childvcf = "{0}/{1}_child{2}.filtered_variants_snpEff_snpSift_Cosmicv76_GoNLv5".format(famdir, family, childnr)
	childvcf_fnameext = "{0}.recode.vcf".format(childvcf)
	if not os.path.isfile(childvcf_fnameext) or args.OVERWRITE == "true":
		command = "vcftools --remove-indv {0} --remove-indv {1} --gzvcf {2} --out {3} --recode --recode-INFO-all".format(sampleids[idstoremove_child[0]], sampleids[idstoremove_child[1]], vcffile, childvcf)
		os.system(command)
		print "Removed samples from a vcf for family: {0}".format(family)
	else:
		print "vcf with removed samples already exists for child{0} in family: {1} and overwrite is false".format(childnr, family)

	#Filter vcf for callable loci
	callablelocibed = "{0}/CallableLoci/Output/{1}_child{2}_intersect_CallableLoci.bed".format(outdir, family, childnr)
	outvcf = "{0}/Callableregion_NrMendelViols/CallableVCF/{1}_filtered_callable_annotated".format(outdir, sample)
	callablevcf = "{0}.recode.vcf".format(outvcf)
	if not os.path.isfile(callablevcf) or args.OVERWRITE == "true":
		command = "vcftools --vcf {0} --bed {1} --out {2} --recode --recode-INFO-all".format(childvcf_fnameext, callablelocibed, outvcf)
		os.system(command)
		print "Created one vcf file with only variants in callable regions"
	else:
		print "vcf file already filtered for 'callable' loci for child{0} in family: {1} and overwrite is false".format(childnr, family)
	
	#Create ped without the other samples
	ped_original = "{0}/{1}/{2}/{2}.ped".format(run, vcfdir, family)
	pattern = "{0}|{1}".format(sampleids[idstoremove_child[0]], sampleids[idstoremove_child[1]])
	ped_child = "{0}/{1}-child{2}-DNA-1.ped".format(famdir, family, childnr)
	if not os.path.isfile(ped_child) or args.OVERWRITE == "true":
		os.system("grep -Ev '{0}' {1} > {2}".format(pattern, ped_original, ped_child))

	#Run phase by transmission
	menviol = "{0}/{1}_child{2}.MendelViol".format(famdir, family, childnr)
	phasedvcf = "{0}/{1}_child{2}.phased.vcf.gz".format(famdir, family, childnr)
	if not os.path.isfile(menviol) or not os.path.isfile(phasedvcf) or args.OVERWRITE == "true":
		command = "java -Xmx2g -jar {0} -T PhaseByTransmission -V {1} -R {2} --MendelianViolationsFile {3} -o {4} -ped {5} --DeNovoPrior {6}".format(jar, callablevcf, genome, menviol, phasedvcf, ped_child, args.MUTATIONPRIOR)
		os.system(command)
		print "Ran PBT"
	else:
		print "Phase by transmission has already been run for child{0} in family: {1} and overwrite is false".format(childnr, family)

	#Remove missing GTs
	Fullgt_menviol = "{0}/{1}_child{2}_fullgenotypes.MendelViol".format(famdir, family, childnr)
	if not os.path.isfile(Fullgt_menviol) or args.OVERWRITE == "true":
		with open(menviol, "r") as MVfile:
			with open(Fullgt_menviol, "w") as outfile:
				header_line = next(MVfile)
				outfile.write(header_line)
				for line in MVfile:
					line = line.rstrip()
					line = line.split("\t")
					mothergt = line[5]
					fathergt = line[9]
					childgt = line[13]
					if "." not in mothergt and "." not in fathergt and "." not in childgt:
						line = "\t".join(line)
						outfile.write("{0}\n".format(line))
		print "Finished removing GTs for child{0} in family: {1}".format(childnr, family)
	else:
		print "Missing GTs have already been removed for child{0} in family: {1} and overwrite is false".format(childnr, family)			
	
	#Intersect mendelviol file with phased vcf
	denovo = "{0}/Callableregion_NrMendelViols/Denovo/{1}_denovo".format(outdir, sample)
	fullname_denovo = "{0}.recode.vcf".format(denovo)
	if not os.path.isfile(fullname_denovo) or args.OVERWRITE == "true":
		command = "vcftools --gzvcf {0} --positions {1} --out {2} --recode --recode-INFO-all".format(phasedvcf, Fullgt_menviol, denovo)
		os.system(command)
		print "Intersected the phased vcf with the mendelviol file for child{0} in family: {1}".format(childnr, family)
	else:
		print "The phased vcf and the mendelviol file have already been intersected for child{0} in family: {1} and overwrite is false".format(childnr, family)

	print "Fully completed run for child{0}".format(childnr)


###End timer
timeend = datetime.datetime.now().strftime("%Y-%m-%d-%H-%M")
print "Finished running script at: {0}".format(timeend)
end = timer()
execution_time = int(end - start)
execution_time = datetime.timedelta(seconds = execution_time)
print "Elapsed time: {0}".format(execution_time)
