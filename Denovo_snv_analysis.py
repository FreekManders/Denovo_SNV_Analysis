import argparse
import os
import sys
import re
import datetime
from timeit import default_timer as timer

#Modules
python = "python/2.7.10"
R = "R/3.4.1"
bedtools = "bed/bedtools/2.25.0"

#Start timer
start = timer()
time = datetime.datetime.now().strftime("%Y-%m-%d-%H-%M")
print "Start running pipeline at: {0}".format(time)


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
	
	
def create_bash(h_rt, h_vmem, scriptname, beas, hold_jid = "", log_script, modules, language, arguments, pipeline_path, output_path, time, log, email):
	"""Creates and runs a bash submission script.
	
	As its input you need to provide the submission variables (h_rt ect.) as well as the name of the script itself, its arguments and the language in which the script was written"""
	
	scriptnamenoextension = os.path.splitext(scriptname)[0]
	bash_fname = "{0}/Jobs/{1}_{2}.sh".format(output_path, scriptnamenoextension, time)
	with open(bash_fname, "w") as bash_file:
		
		#Create the general settings of the bash script
		bash_file.write("#!/bin/bash\n#$ -S /bin/bash\n#$ -l h_rt={0}\n#$ -l h_vmem={1}\n#$ -cwd\n#$ -N {2}\n".format(h_rt, h_vmem, scriptnamenoextension))
		
		if log_script.lower() == "true":
			bash_file.write("#$ -o {0}/Logs/{1}_{2}_output\n#$ -e {0}/Logs/{1}_{2}_errors\n".format(output_path, scriptnamenoextension, time))
		else:
			bash_file.write("#$ -o /dev/null\n#$ -e /dev/null\n")
		
		if hold_jid is not "":
			bash_file.write("#$ -hold_jid {0}\n".format(hold_jid))
		
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


#Get location of the pipeline
pathname = os.path.dirname(sys.argv[0])        
pipeline_path = os.path.abspath(pathname)

#Parse arguments
parser = argparse.ArgumentParser(description = "This script runs the entire autosomal denovo snv analyis pipeline. Output from the UMCU IAP pipeline is used as its input.")
parser.add_argument("-i", "--ini", default = "/hpc/cog_bioinf/cuppen/project_data/Freek_SNV/pipeline/Denovo_snv_analysis.ini", help = "The input ini file")
args = parser.parse_args()

#Parse the ini file
ini_dict = ini_parser(args.ini)
output_path = ini_dict["OUTPUT_PATH"]
filelist = ini_dict["FILELIST"]
log = ini_dict["OVERVIEW_LOG"].lower()
email = ini_dict["EMAIL"]
overwrite = ini_dict["OVERWRITE"].lower()
	
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

#Run the different parts of the pipeline
if ini_dict["WGS_QC"].lower() == "true":
	arguments = "--FILELIST {0} --OUTPUT_PATH {1} --OVERWRITE {2}".format(filelist, output_path, overwrite)
	create_bash(h_rt = ini_dict["WGS_QC_TIME"], h_vmem = ini_dict["WGS_QC_MEM"], scriptname = "WGS_QC_plots.R", beas = ini_dict["WGS_QC_BEAS"], log_script = ini_dict["WGS_QC_LOG"], modules = [R], language = "Rscript", arguments = arguments, pipeline_path = pipeline_path, output_path = output_path, time = time, log = log, email = email)

if ini_dict["CallableLoci"].lower() == "true":
	arguments = "--FILELIST {0} --OUTPUT_PATH {1} --OVERWRITE {2}".format(filelist, output_path, overwrite)
	create_bash(h_rt = ini_dict["CallableLoci_TIME"], h_vmem = ini_dict["CallableLoci_MEM"], scriptname = "GetCallableLoci.py", beas = ini_dict["CallableLoci_BEAS"], log_script = ini_dict["CallableLoci_LOG"], modules = [python, bedtools], language = "python", arguments = arguments, pipeline_path = pipeline_path, output_path = output_path, time = time, log = log, email = email)

#End timer
timeend = datetime.datetime.now().strftime("%Y-%m-%d-%H-%M")
print "Finished submitting pipeline at: {0}".format(timeend)
end = timer()
execution_time = int(end - start)
execution_time = datetime.timedelta(seconds = execution_time)
print "Elapsed time: {0}".format(execution_time)

		
