import argparse
import os
import sys
import re
import datetime
from timeit import default_timer as timer
import pandas as pd

####____________________Start timer______________________####
start = timer()
time = datetime.datetime.now().strftime("%Y-%m-%d-%H-%M")
print "Start running pipeline at: {0}".format(time)

####____________________Function to parse ini files_______####
def ini_parser(ini_fname):
	r"""Parses a ini file with a key\tvalue format and returns a dictionary"""
	if not os.path.isfile(ini_fname): raise IOError("The ini file could not be read. This has caused an error.")
	ini_dict = {}
	with open(ini_fname) as ini:
		for line in ini:
			line = line.strip()
			if line.startswith("#") or line == "":
				continue
			if "\t" in line:
				key, value = line.split("\t")
			else:
				key = line
				value = ""
			ini_dict[key] = value
	return ini_dict

####_________Annotation functions_______####
def annotategnomAD(filelist, output_path, overwrite, snpsift, gnomad_folder, time):
	
	folder = "{0}/Callableregion_NrMendelViols/TrueDenovo/SNVs".format(output_path, gnomad_folder)
	gnomADfolder = gnomad_folder
	gnomADhdr = "{0}/gnomADannotation.hdr".format(gnomADfolder)
	annotatedf = "{0}/Tables/gnomADannotated.txt".format(folder)
	jobfolder = "{0}/Jobs/gnomADjobs".format(output_path)
	
	###Read the filelist
	dffiles = pd.read_csv(filelist, sep = "\t", header = 0)
	childrowsf = dffiles["Name"].str.contains("Child|child")
	dffiles = dffiles[childrowsf]
	dffiles["vcf"] = folder + "/" + dffiles["Sample"] + "_Truedenovo.vcf"
	dffiles["vcftempout"] = folder + "/Temp/" + dffiles["Sample"] + "_Truedenovo.vcf"

	###Bash settings
	bashsettings = "#!/bin/bash\n#$ -S /bin/bash\n#$ -l h_rt=00:10:00\n#$ -l h_vmem=1G\n#$ -cwd\n#$ -N annognomAD\n#$ -hold_jid annotateDANN\n#$ -o {0}/Logs/annotategnomAD_{1}_output\n#$ -e {0}/Logs/annotategnomAD_{1}_errors\n#$ -M fmanders@umcutrecht.nl\n#$ -m as\nmodule load vcfbcf/bcftools/1.3\n".format(output_path, time)

	bashsettingsjoin = "#!/bin/bash\n#$ -S /bin/bash\n#$ -l h_rt=00:05:00\n#$ -l h_vmem=7G\n#$ -cwd\n#$ -N gnomADAnnotationJoining\n#$ -hold_jid annognomAD\n#$ -o {0}/Logs/annotategnomAD_{1}_output\n#$ -e {0}/Logs/annotategnomAD_{1}_errors\n#$ -M fmanders@umcutrecht.nl\n#$ -m as\n".format(output_path, time)


	#Remove samples that have already been annotated from to do list
	if not overwrite == "true" and os.path.isfile(annotatedf):
		annotated = pd.read_csv(annotatedf, header = None)
		alreadyrun = dffiles["Sample"].isin(annotated[0])
		dffiles = dffiles[~alreadyrun]
		nr2run = dffiles.shape[0]

	###Put filepaths in a nested list. This list conains a list with filepaths for each sample.
	lfiles = dffiles[["Sample", "vcf", "vcftempout"]].values.tolist()


	###Function to annotate samples with gnomAD frequencies
	def AnnotategnomAD_sample(lfiles_sample):
		sample, vcf, vcftempout = lfiles_sample
		if os.path.exists(vcf):
			command = "java -Xmx4G -jar {0} split {1}".format(snpsift, vcf)
			os.system(command)
		
			#Create a bashfile per sample
			bashname = "{0}/Annotate_gnomAD_{1}_{2}.sh".format(jobfolder, sample, time)
			with open(bashname, "w") as bash:
				bash.write(bashsettings)
		
			#Create a command per chromosome
			for i in range(1,23):
				vcfchrom = vcf.rstrip(".vcf")
				vcfchrom = "{0}.{1}.vcf".format(vcfchrom, i)
				vcfchromtempout = "{0}/Temp/{1}_{2}_temp.vcf".format(folder, sample, i)
				gnomadchrom = "{0}/gnomADpass_snv_table_{1}.txt.gz".format(gnomADfolder, i)

			
				if os.path.exists(vcfchrom):
					command1 = "bcftools annotate -a {0} -c CHROM,POS,REF,ALT,gnomAD_AC,gnomAD_AF,gnomAD_AN {1} -h {2} -o {3}".format(gnomadchrom, vcfchrom, gnomADhdr, vcfchromtempout)
					command2= "rm {0}".format(vcfchrom)
					fullbash = "{0}\n{1}\n".format(command1, command2)
					with open(bashname, "a") as bash:
						bash.write(fullbash)
			#Finish making the bashfile
			with open(bashname, "a") as bash:
				bash.write("module unload vcfbcf/bcftools/1.3")
			os.system("qsub {0}".format(bashname))
		
			#Write bashfile for joining annotated vcfs
			bashname = "{0}/Annotate_gnomAD_joining_{1}_{2}.sh".format(jobfolder, sample, time)
			command1 = """vcfchromanno="$(ls {0}/Temp/{1}_*_temp.vcf | sort -V)" """.format(folder, sample)
			command2 = """java -Xmx4G -jar {0} split -j  ${{vcfchromanno}} > {1}""".format(snpsift, vcf)
			command3 = """rm ${vcfchromanno}"""
			fullbash = "{0}\n{1}\n{2}\n{3}".format(bashsettingsjoin, command1, command2, command3)
			with open(bashname, "w") as bash:
						bash.write(fullbash)
			os.system("qsub {0}".format(bashname))

			#Add sample to list of annotated samples
			with open(annotatedf, "a") as annotated:
				annotated.write("{0}\n".format(sample))
			
	###Run annotation function
	map(AnnotategnomAD_sample, lfiles)

	###Removes duplicates from the list of already annotated samples.
	with open(annotatedf, "r") as annotated:
		uniqlines = sorted(set(annotated.readlines()))

	with open(annotatedf, "w") as annotated:
		annotated.writelines(uniqlines)

def annotateHMFPON(filelist, output_path, overwrite, snpsift, time, hmfpon_folder):
	#Create the file paths.
	filelist = "/hpc/cog_bioinf/cuppen/project_data/Freek_SNV/FilelistFreek.txt"
	folder = "{0}/Callableregion_NrMendelViols/TrueDenovo/SNVs".format(output_path)
	jobfolder = "{0}/Jobs/HMFPONjobs".format(output_path)
	HMFPON = "{0}/PON_v2.0_pass_snv.vcf.gz".format(hmfpon_folder)
	annotatedf = "{0}/Callableregion_NrMendelViols/TrueDenovo/SNVs/Tables/HMF-PONannotated.txt".format(output_path)

	#Read the file list
	dffiles = pd.read_csv(filelist, sep = "\t", header = 0)
	childrowsf = dffiles["Name"].str.contains("Child|child")
	dffiles = dffiles[childrowsf]
	dffiles["vcf"] = folder + "/" + dffiles["Sample"] + "_Truedenovo.vcf"
	dffiles["vcftempout"] = folder + "/Temp/" + dffiles["Sample"] + "_Truedenovo.vcf"

	bashsettings = "#!/bin/bash\n#$ -S /bin/bash\n#$ -l h_rt=00:30:00\n#$ -l h_vmem=10G\n#$ -cwd\n#$ -N annoHMFPON\n#$ -hold_jid gnomADAnnotationJoining\n#$ -o {0}/Logs/annotateHMF-PON_{1}_output\n#$ -e {0}/Logs/annotateHMF-PON_{1}_errors\n#$ -M fmanders@umcutrecht.nl\n#$ -m as\n".format(output_path, time)

	#Remove samples that have already been annotated from to do list
	if not overwrite == "true" and os.path.isfile(annotatedf):
		annotated = pd.read_csv(annotatedf, header = None)
		alreadyrun = dffiles["Sample"].isin(annotated[0])
		dffiles = dffiles[~alreadyrun]


	#Put filepaths in a nested list. This list conains a list with filepaths for each sample.
	lfiles = dffiles[["Sample", "vcf", "vcftempout"]].values.tolist()

	#Annotate the samples
	def AnnotateHMFPON_sample(lfiles_sample):
		sample, vcf, vcftempout = lfiles_sample
		if os.path.exists(vcf):
			bashname = "{0}/Annotate_HMFPON_{1}_{2}.sh".format(jobfolder, sample, time)
			command1 = "java -Xmx7G -jar {0} annotate -info PON_COUNT {1} {2} > {3}".format(snpsift, HMFPON, vcf, vcftempout)
			command2 = "mv {0} {1}".format(vcftempout, vcf)
			fullbash = "{0}\n{1}\n{2}".format(bashsettings, command1, command2)
			with open(bashname, "w") as bash:
				bash.write(fullbash)
			os.system("qsub {0}".format(bashname))
		
			#Add sample to list of annotated samples
			with open(annotatedf, "a") as annotated:
				annotated.write("{0}\n".format(sample))

	map(AnnotateHMFPON_sample, lfiles)

	#Removes duplicates from the list of already annotated samples.
	with open(annotatedf, "r") as annotated:
		uniqlines = sorted(set(annotated.readlines()))

	with open(annotatedf, "w") as annotated:
		annotated.writelines(uniqlines)

def Annotatevcfwithbed(filelist, output_path, overwrite, snpsift, time, bed_folders, id2cell_chromhmm):
	
	#Create some filepaths based on the arguments.
	folder = "{0}/Callableregion_NrMendelViols/TrueDenovo/SNVs".format(output_path)
	jobfolder = "{0}/Jobs/bedjobs".format(output_path)
	annotatedf = "{0}/Tables/bedannotated.txt".format(folder)
	bed_folders = bed_folders.split(",")
	
	#Read file list
	dffiles = pd.read_csv(filelist, sep = "\t", header = 0)
	childrowsf = dffiles["Name"].str.contains("Child|child")
	dffiles = dffiles[childrowsf]
	dffiles["vcf"] = folder + "/" + dffiles["Sample"] + "_Truedenovo.vcf"
	dffiles["vcftempout"] = folder + "/Temp/" + dffiles["Sample"] + "_Truedenovo.vcf"
	dffiles["vcftempout2"] = folder + "/Temp/" + dffiles["Sample"] + "_Truedenovo2.vcf"

	bashsettings = "#!/bin/bash\n#$ -S /bin/bash\n#$ -l h_rt=00:40:00\n#$ -l h_vmem=10G\n#$ -cwd\n#$ -N annobed\n#$ -hold_jid annoHMFPON\n#$ -o {0}/Logs/annotatebeds_{1}_output\n#$ -e {0}/Logs/annotatebeds_{1}_errors\n#$ -M fmanders@umcutrecht.nl\n#$ -m as\nmodule load vcfbcf/bcftools/1.3\n".format(output_path, time)

	#Find bedfiles to annotate.
	beds = []
	for bedfolder in bed_folders:
		for file in os.listdir(bedfolder):
			if file.endswith(".bed.gz"):
				bedf = os.path.join(bedfolder, file)
				beds.append(bedf)

	#Create annotation names based on the first part of the bedfile name (and a key2name file(For the ChromHMM bed files).
	id2cell = pd.read_csv(id2cell_chromhmm, sep = "\t", header = 0)

	annots = []
	for bed in beds:
		annot = bed.split("/")[-1]
		annot = re.sub(r"(\.|_).*", "", annot)
		row = id2cell.loc[id2cell["id"] == annot]
		if row.shape[0] == 1:
			annot = row.iloc[0]["celltype"]
		annots.append(annot)

	#Create the header filenames to use with the bedfiles for annotation
	hdrfs = []
	for i in range(0, len(beds)):
		bed = beds[i]
		hdrf = re.sub(".bed.gz", ".hdr", bed)
		hdrfs.append(hdrf)

	#Write the actual header files to use with the bedfiles for annotation
	for i in range(0, len(beds)):
		hdrf = hdrfs[i]
		annot = annots[i]
		line = """##INFO=<ID={0},Number=1,Type=String,Description="Shows the function of the variant in {0}">""".format(annot)
		with open(hdrf, "w") as hdr:
			hdr.write(line)

	#Remove samples that have already been annotated from to do list
	if not overwrite == "true" and os.path.isfile(annotatedf):
		annotated = pd.read_csv(annotatedf, header = None)
		alreadyrun = dffiles["Sample"].isin(annotated[0])
		dffiles = dffiles[~alreadyrun]


	#Annotate samples with bed files
	for index, row in dffiles.iterrows():
		sample, vcf, vcftempout, vcftempout2 = row[["Sample", "vcf", "vcftempout", "vcftempout2"]].values.tolist()
		if os.path.exists(vcf):
		
			#Write start of bash
			bashname = "{0}/{1}_{2}.sh".format(jobfolder, sample, time)
			commandstart = "cp {0} {1}\n".format(vcf, vcftempout)
			with open(bashname, "w") as bash:
				bash.write(bashsettings)
				bash.write(commandstart)
		
			#For each bedfile	
			for i in range(0, len(beds)):
				bed = beds[i]
				annot = annots[i]
				hdrf = hdrfs[i]
				command1 = "bcftools annotate -a {0} -c CHROM,FROM,TO,{1} {2} -h {3} -o {4}\n".format(bed, annot, vcftempout, hdrf, vcftempout2)
				command2 = "mv {0} {1}\n".format(vcftempout2, vcftempout)
				with open(bashname, "a") as bash:
					bash.write(command1)
					bash.write(command2)
		
			#Write end of bash
			commandend = "mv {0} {1}\n".format(vcftempout, vcf)
			with open(bashname, "a") as bash:
				bash.write(commandend)
				bash.write("module unload vcfbcf/bcftools/1.3")
			os.system("qsub {0}".format(bashname))
		
			#Add sample to list of annotated samples
			with open(annotatedf, "a") as annotated:
				annotated.write("{0}\n".format(sample))
		

	
	#Removes duplicates from the list of already annotated samples.
	with open(annotatedf, "r") as annotated:
		uniqlines = sorted(set(annotated.readlines()))

	with open(annotatedf, "w") as annotated:
		annotated.writelines(uniqlines)

####______________________Get location of the pipeline____####
pathname = os.path.dirname(sys.argv[0])        
pipeline_path = os.path.abspath(pathname)

####______________________Parse arguments_________________####
parser = argparse.ArgumentParser(description = "This script runs the entire autosomal denovo snv analyis pipeline. Output from the UMCU IAP pipeline is used as its input.")
parser.add_argument("-i", "--ini", default = "/hpc/cog_bioinf/cuppen/project_data/Freek_SNV/pipeline/Denovo_snv_analysis.ini", help = "The input ini file")
args = parser.parse_args()

####_______________________Parse the ini file_______________####
ini_dict = ini_parser(args.ini)

#Main configurations
output_path = ini_dict["OUTPUT_PATH"]
filelist = ini_dict["FILELIST"]
log = ini_dict["OVERVIEW_LOG"].lower()
email = ini_dict["EMAIL"]
overwrite = ini_dict["OVERWRITE"].lower()

#Modules. These can be loaded via the module system.
python = ini_dict["PYTHON"]
R = ini_dict["R"]
bedtools = ini_dict["BEDTOOLS"]
vcftools = ini_dict["VCFTOOLS"]
tabix = ini_dict["TABIX"]
bcftools = ini_dict["BCFTOOLS"]

#Other software. These are often .jar fileS.
gatk = ini_dict["GATK"]
snpsift = ini_dict["SNPSIFT"]
denovogear = ini_dict["DENOVOGEAR"]

#Data. This is data, that wasn't generated by the IAP or this pipeline.
genome = ini_dict["GENOME"]

####__________________General preparations____________________####
#Create the directories needed in the output_path.
with open(os.path.join(pipeline_path, "dirlist.txt")) as dirlist:
	for line in dirlist:
		dir_path = line.strip()
		dir_path = os.path.join(output_path, dir_path)
		if not os.path.isdir(dir_path):
			os.mkdir(dir_path)

#initialize a log file of all the modules that have been used.
if log == "true":
	modulesusedlog_fname = "{0}/Logs/modules_used_{1}.txt".format(output_path, time)
	with open(modulesusedlog_fname, "w") as f:
		pass

#Standard arguments for the scripts.
std_arguments = "--FILELIST {0} --OUTPUT_PATH {1}".format(filelist, output_path)

####_________________________Function to create bash files.__________________________________####
def create_bash(mode, scriptname, modules, language, arguments = std_arguments, pipeline_path = pipeline_path, output_path = output_path, time = time, log = log, email = email, hold_jid = "", threads = 1, nr = ""):
	"""Creates and runs a bash submission script.
	
	As its input you need to provide the mode to run as well as the name of the script itself, its arguments and the language in which the script was written. The function also uses several global variables as default values."""
	
	#Look up some of the ini settings
	if ini_dict[mode].lower() != "true":
		return
		
	h_rt = ini_dict["{0}_TIME".format(mode)]
	h_vmem = ini_dict["{0}_MEM".format(mode)]
	beas = ini_dict["{0}_BEAS".format(mode)]
	log_script = ini_dict["{0}_LOG".format(mode)]
	
	#start writing the bash
	scriptnamenoextension = os.path.splitext(scriptname)[0]
	bash_fname = "{0}/Jobs/{1}_{2}{3}.sh".format(output_path, scriptnamenoextension, time, nr)
	with open(bash_fname, "w") as bash_file:
		
		#Create the general settings of the bash script
		bash_file.write("#!/bin/bash\n#$ -S /bin/bash\n#$ -l h_rt={0}\n#$ -l h_vmem={1}\n#$ -cwd\n#$ -N {2}\n".format(h_rt, h_vmem, scriptnamenoextension))
		
		if log_script.lower() == "true":
			bash_file.write("#$ -o {0}/Logs/{1}_{2}_output\n#$ -e {0}/Logs/{1}_{2}_errors\n".format(output_path, scriptnamenoextension, time))
		else:
			bash_file.write("#$ -o /dev/null\n#$ -e /dev/null\n")
		
		if hold_jid is not "":
			bash_file.write("#$ -hold_jid {0}\n".format(hold_jid))
		
		if threads != 1:
			bash_file.write("#$ -pe threaded {0}\n".format(threads))
		
		if email is not "":
			bash_file.write("#$ -M {0}\n".format(email))
		
		if beas is not "":
			bash_file.write("#$ -m {0}\n".format(beas))
		
		#Log the time the scripts starts running
		if log == "true":
			bash_file.write("""echo "running {0} on $(date)" >> {1}\n""".format(scriptnamenoextension, modulesusedlog_fname))
	
		#module loading
		for module in modules:
			bash_file.write("module load {0}\n".format(module))

		#Create the call to the script that needs to be run
		script_abspath = os.path.join(pipeline_path, scriptname)
		bash_file.write("{0} {1} {2}\n".format(language, script_abspath, arguments))
	
		#module unloading
		for module in modules:
			bash_file.write("module unload {0}\n".format(module))
	
		#Log the time the script finished running.
		if log == "true":
			bash_file.write("""echo "Finished running {0} on $(date)" >> {1}\n""".format(scriptnamenoextension, modulesusedlog_fname))
		
	#Run the bash file
	os.system("qsub {0}".format(bash_fname))


####__________Run the different parts of the pipeline._______####	

#WGS_QC
create_bash(mode = "WGS_QC", scriptname = "WGS_QC_plots.R", modules = [R], language = "Rscript")

#CallableLoci
arguments = "{0} --OVERWRITE {1}".format(std_arguments, overwrite)
create_bash(mode = "CallableLoci", scriptname = "GetCallableLoci.py", modules = [python, bedtools], language = "python", arguments = arguments)

#Callableregions_and_runPbt
arguments = "{0} --GATK {1} --GENOME {2} --THREADS {3} --OVERWRITE {4}".format(std_arguments, gatk, genome, ini_dict["Callableregions_and_runPbt_THREADS"], overwrite)
create_bash(mode = "Callableregions_and_runPbt", scriptname = "Callableregions_and_runPbt.py", threads = ini_dict["Callableregions_and_runPbt_THREADS"], hold_jid = "GetCallableLoci", modules = [python, vcftools], language = "python", arguments = arguments)

#Remove_missingGTs_MendelViol
arguments = "{0} --OVERWRITE {1}".format(std_arguments, overwrite)
create_bash(mode = "Remove_missingGTs_MendelViol", scriptname = "Remove_missingGTs_MendelViol.py", hold_jid = "Callableregions_and_runPbt", modules = [python], language = "python", arguments = arguments)

#Intersect_mendelviol_phasedvcf
arguments = "{0} --OVERWRITE {1}".format(std_arguments, overwrite)
create_bash(mode = "Intersect_mendelviol_phasedvcf", scriptname = "Intersect_mendelviol_phasedvcf.py", hold_jid = "Remove_missingGTs_MendelViol", modules = [python, vcftools], language = "python", arguments = arguments)

#Callableregions_and_runPbt_3kids
families_3kids = ini_dict["Callableregions_and_runPbt_3kids_FAMILY"]
families_3kids = families_3kids.split(",")
for i in range(len(families_3kids)):
	family = families_3kids[i]
	nr = "_nr{0}".format(i+1)
	arguments = "{0} --GATK {1} --GENOME {2} --FAMILY {3} --OVERWRITE {4}".format(std_arguments, gatk, genome, family, overwrite)
	create_bash(mode = "Callableregions_and_runPbt_3kids", scriptname = "Callableregions_and_runPbt_3kids.py", hold_jid = "GetCallableLoci", modules = [python, vcftools], language = "python", arguments = arguments, nr = nr)

#FilterDenovo
arguments = "{0} --SNPSIFT {1} --OVERWRITE {2}".format(std_arguments, snpsift, overwrite)
create_bash(mode = "FilterDenovo", scriptname = "FilterDenovo.py", hold_jid = "Intersect_mendelviol_phasedvcf,Callableregions_and_runPbt_3kids", modules = [python, vcftools], language = "python", arguments = arguments)

#denovogear
arguments = "{0} --SNPSIFT {1} --DENOVOGEAR {2} --OVERWRITE {3}".format(std_arguments, snpsift, denovogear, overwrite)
create_bash(mode = "denovogear", scriptname = "denovogear.py", hold_jid = "FilterDenovo", modules = [python, bcftools, tabix], language = "python", arguments = arguments)

#manualcheckdenovos
arguments = "{0} --SNPSIFT {1} --GATK {2} --GENOME {3}".format(std_arguments, snpsift, gatk, genome)
create_bash(mode = "manualcheckdenovos", scriptname = "manualcheckdenovos.py", hold_jid = "denovogear", modules = [python, vcftools], language = "python", arguments = arguments)

#Filtermanualfps
create_bash(mode = "Filtermanualfps", scriptname = "Filtermanualfps.py", hold_jid = "manualcheckdenovos", modules = [python, vcftools], language = "python")

#annotateDANN
arguments = "{0} --OVERWRITE {1} --DANN_FOLDER {2} --THREADS {3}".format(std_arguments, overwrite, ini_dict["DANN_FOLDER"], ini_dict["annotateDANN_THREADS"])
create_bash(mode = "annotateDANN", scriptname = "annotateDANN.py", hold_jid = "Filtermanualfps", modules = [python, bcftools], language = "python", arguments = arguments, threads = ini_dict["annotateDANN_THREADS"])

#annotategnomAD
if ini_dict["annotategnomAD"].lower() == "true":
	annotategnomAD(filelist = filelist, output_path = output_path, overwrite = overwrite, gnomad_folder = ini_dict["GNOMAD_FOLDER"], snpsift = snpsift, time = time)

#annotateHMFPON
if ini_dict["annotateHMFPON"].lower() == "true":
	annotateHMFPON(filelist = filelist, output_path = output_path, overwrite = overwrite, hmfpon_folder = ini_dict["HMFPON_FOLDER"], snpsift = snpsift, time = time)

#annotate with bed files
if ini_dict["annotatevcfwithbed"].lower() == "true":
	Annotatevcfwithbed(filelist = filelist, output_path = output_path, overwrite = overwrite, snpsift = snpsift, time = time, bed_folders = ini_dict["ANNOTATION_BED_FOLDERS"], id2cell_chromhmm = ini_dict["ID_2_CELLTYPE_CHROMHMM"])

#Filterparentalorigin
arguments = "{0} --SNPSIFT {1} --OVERWRITE {2}".format(std_arguments, snpsift, overwrite)
create_bash(mode = "Filterparentalorigin", scriptname = "Filterparentalorigin.py", hold_jid = "annobed", modules = [python, vcftools], language = "python", arguments = arguments)

####____________________End timer_____________________####
timeend = datetime.datetime.now().strftime("%Y-%m-%d-%H-%M")
print "Finished submitting pipeline at: {0}".format(timeend)
end = timer()
execution_time = int(end - start)
execution_time = datetime.timedelta(seconds = execution_time)
print "Elapsed time: {0}".format(execution_time)

		
