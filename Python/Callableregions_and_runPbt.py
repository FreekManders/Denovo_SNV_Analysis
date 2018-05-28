from multiprocessing import Pool as ThreadPool
import os
import argparse
import datetime
from timeit import default_timer as timer
import subprocess
import signal

###Start timer
start = timer()
time = datetime.datetime.now().strftime("%Y-%m-%d-%H-%M")
print "Start running script at: {0}".format(time)

###Parse arguments
parser = argparse.ArgumentParser(description = "This script filters vcf files on 'callable' loci and then runs it through PhaseByTransmission.")
parser.add_argument("-f", "--FILELIST", default = "/hpc/cog_bioinf/cuppen/project_data/Freek_SNV/FilelistFreek.txt", help = "The file list")
parser.add_argument("-o", "--OUTPUT_PATH", default = "/hpc/cog_bioinf/cuppen/project_data/Freek_SNV", help = "The output path")
parser.add_argument("-w", "--OVERWRITE", default = "False", help = "Overwrite the results if they already exist.")
parser.add_argument("-g", "--GENOME", default = "/hpc/cog_bioinf/GENOMES/Homo_sapiens.GRCh37.GATK.illumina/Homo_sapiens.GRCh37.GATK.illumina.fasta", help = "The used genome in fasta format.")
parser.add_argument("--GATK", default = "/hpc/cog_bioinf/common_scripts/GenomeAnalysisTK-3.4-46/GenomeAnalysisTK.jar", help = "The GATK jar.")
parser.add_argument("--MUTATIONPRIOR", default = "1.0E-4", help = "The mutation prior that is being used.")
parser.add_argument("-t", "--THREADS", default = 1, type = int, help = "Number of threads to use")
args = parser.parse_args()

filelist = args.FILELIST
genome = args.GENOME
jar = args.GATK
outdir = "{0}/Callableregion_NrMendelViols".format(args.OUTPUT_PATH)
callablelocidir = "{0}/CallableLoci/Output".format(args.OUTPUT_PATH)

###Main body of script
commands = []
families = []
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
			run = columns[1]
			
			###Find the correct vcf file
			if os.path.isdir("{0}/vcf_new_filter".format(run)):
				vcfdir = "vcf_new_filter"
			else:
				vcfdir = "vcf"
			vcffile = "{0}/{1}/{2}/{2}.filtered_variants_snpEff_snpSift_Cosmicv76_GoNLv5.vcf.gz".format(run, vcfdir, family)
			if not os.path.isfile(vcffile):
				vcffile = "{0}/{1}/{2}/{2}.filtered_variants_snpEff_snpSift_Cosmicv76_GoNLv5.vcf".format(run, vcfdir, family)
				if not os.path.isfile(vcffile):
					print "vcf file does not exists."
					continue
			
			sampleids = subprocess.check_output(["vcf-query", "-l", vcffile], preexec_fn=lambda:signal.signal(signal.SIGPIPE, signal.SIG_DFL))
			sampleids = sampleids.split("\n")
			del sampleids[-1] #Remove empty newline
			if len(sampleids) != 3:
				continue
			
			###Create command to filter vcf files based on the CallableLoci.bed files.
			callablelocibed = "{0}/{1}_intersect_CallableLoci.bed".format(callablelocidir, family)
			outvcf = "{0}/CallableVCF/{1}_filtered_callable_annotated".format(outdir, sample)
			callablevcf = "{0}.recode.vcf".format(outvcf)
			if os.path.isfile(callableLocibed) and (not os.path.isfile(callablevcf) or args.OVERWRITE == "true"):
				command1 = "vcftools --gzvcf {0} --bed {1} --out {2} --recode --recode-INFO-all; echo Finished filtering the vcf file for family: {3};".format(vcffile, callablelocibed, outvcf, family)
				print "Created command to make vcf file with only variants in callable regions for family: {0}".format(family)
			else:
				print "vcf file for family: {0} already created and overwrite is false".format(family)
				command1 = ""
			
			###Create command to run PhaseByTransmission on the filtered vcf files.
			ped = "{0}/{1}/{2}/{2}.ped".format(run, vcfdir, family)
			menviol = "{0}/PBT/{1}.MendelViol".format(outdir, family)
			phasedvcf = "{0}/PBT/{1}.phased.vcf.gz".format(outdir, family)
			if not os.path.isfile(menviol) or not os.path.isfile(phasedvcf) or args.OVERWRITE == "true":
				command2 = "java -Xmx10g -jar {0} -T PhaseByTransmission -V {1} -R {2} --MendelianViolationsFile {3} -o {4} -ped {5} --DeNovoPrior {6}; echo Finished with PhaseByTransmission for family: {7}".format(jar, callablevcf, genome, menviol, phasedvcf, ped, args.MUTATIONPRIOR, family)
				print "Created command to run PhaseByTransmission for family {0}".format(family)
			else:
				print "Phase by transmission has already been run for family: {0} and overwrite is false".format(family)
				command2 = "" 
			combicommand = "{0} {1}".format(command1, command2)
			commands.append(combicommand)

###Run the commands on multiple threads
pool = ThreadPool(args.THREADS)
results = pool.map(os.system, commands)
pool.close()
pool.join()

###End timer
timeend = datetime.datetime.now().strftime("%Y-%m-%d-%H-%M")
print "Finished running script at: {0}".format(timeend)
end = timer()
execution_time = int(end - start)
execution_time = datetime.timedelta(seconds = execution_time)
print "Elapsed time: {0}".format(execution_time)


