#######################################################################
# Script: Run FJD analysis for samples stored in basespace or locally
# Author: Lorena de la Fuente 
# Date: 24-01-2019
#######################################################################

holaaaaaaaaaaaaaaaaaaaaaaaaaaaaa
#!/usr/bin/env python

import sys
import os
from glob import glob
import subprocess
import argparse
from collections import Counter
import datetime
import time

now = (datetime.datetime.now()).strftime('%Y_%m_%d_%H_%M_%S')
print("DATE:"+now)
utilitiesPath =  os.path.dirname(os.path.realpath(__file__))+"/utilities/" 
sys.path.insert(0, utilitiesPath)

from concatenation_1 import joinFastq



def sbatch(job_name, folder_out, command, time=4, mem=5, threads=5, mail=None, dep=''):

	if dep != '':
		dep = '--dependency=afterok:{} --kill-on-invalid-dep=yes '.format(dep)

	if mail!=None:
		mailc = "--mail-user={} --mail-type=FAIL".format(mail)
	else:
		mailc = ''

	sbatch_command = "sbatch -J {} -o {}/{}.out -e {}/{}.err {} -t {}:00:00 --account=bioinfofjd_serv --partition=bioinfofjd --mem-per-cpu={}gb --cpus-per-task={} {} {}".format(job_name, folder_out, job_name, folder_out, job_name, mailc, time, mem, threads, dep,command)
	print(sbatch_command)

	sbatch_response = subprocess.check_output(sbatch_command, shell=True)
	job_id = sbatch_response.split(' ')[-1].strip()
	return job_id
			




def main():

	#arguments
	parser = argparse.ArgumentParser(description="FJD genotyping analysis for samples from Basespace")
	parser.add_argument('-i', '--input', help='\t\t Local directory with fastq.gz/bam input files or project name within Basespace', required=True)
	parser.add_argument('-o', '--output', dest="output",help='\t\tOutput directory name.', required=True)
	parser.add_argument('-s', '--samples' , help='\t\tCSV file listing sample names within second column', required=False, default="all")
	parser.add_argument('-n', '--name', help='\t\tName for job. RECOMMENDED', required=False)
	parser.add_argument('-a', '--analysis', help='\t\tType of analysis to run', required=False, choices={"mapping", "snp", "cnv", "all"}, default="snp")
	parser.add_argument('-c', '--cvcf', help='\t\tCombined genotyping. Number of samples must be higher than 2', required=False, action='store_true')
	parser.add_argument('-p', '--panel', help='\t\tBed with panel regions', required=False)
	parser.add_argument('-b', '--basespace', help='\t\tTake samples from Basespace', required=False, action='store_true')
	parser.add_argument('-l', '--local', help='\t\tRun in local', required=False, action='store_true')
	parser.add_argument('-t', '--threads', help='\t\tNumber of threads', type=int, required=False, default=1)
	parser.add_argument('-m', '--mail', help='\t\tMail account', required=False)
	parser.add_argument('-k', '--skipMapping', help='\t\tOption to skip mapping. Input folder must be fastq folder and Output folder the global output folder', action='store_true')
	parser.add_argument('-g', '--genome', help='\t\tLocal directory for genome reference file', required=False)
	parser.add_argument('-e', '--pedigree', help='\t\tPedigree file. Combined genotyping will be run automatically', required=False)
	parser.add_argument('-v', '--version', help="Display program version number.", action='version', version='1.0')
	parser.add_argument('-P', '--pathology', help="Disease group name: resp, digest, skin, musc, gu, pregnancy, perinatal, cong_mal, clinic, infectious, other, neoplasm, blood, endoc, mental, CNS, eye, ear, circ,healthy.", default='healthy')
	parser.add_argument('-I', '--intervals', help="Specified if genomic intervals (panel) over which to operate during SNV calling are used.", action='store_true')
	#parser.add_argument('-C', '--cnv-method', help="Method to call copy number variants (CNVs). Comma-separated list of methods. Methods: ED, CN, C2, PM.", required=False, default="ED")

	args = parser.parse_args()



	# label for current lab

	if args.basespace:
		run = args.input+"_"+now
	else:
		if args.name!=None:
			run = args.name+"_"+now
		else:
			run = now




	## defining stdout and stderr file names 

	stdoutAll="FJD_"+run+".out"
	sys.stdout = open(stdoutAll, 'w', 0)
	stderrAll="FJD_"+run+".err"
	sys.stderr = open(stderrAll, 'w', 0)





	## pipeline starts

	sys.stdout.write("\nFJD ANALYSIS\n")
	sys.stdout.write("\nChecking arguments...\n")






	# checking output dir

	if not os.path.isdir(args.output): 
		sys.stderr.write("ERROR: Output folder '%s' does not exist\n" %(args.output))
		sys.exit()
	else:
		args.output=os.path.realpath(args.output)







	# checking if panel file exists

	if args.panel!= None:
		if not os.path.isfile(args.panel): 
			sys.stderr.write("ERROR: Panel file '%s' does not exist\n" %(args.panel))
			sys.exit()
		else:
			args.panel=os.path.realpath(args.panel)

	else:
		args.panel="genome"
		if args.intervals == True:
			sys.stderr.write("ERROR: Intervals specified but bed file not found.\n")

	

	





	# checking if specified bundle exists

	if args.genome!=None:
		if not os.path.isdir(args.genome): 
			sys.stderr.write("ERROR: Genome bundle folder '%s' does not exist\n" %(args.genome))
			sys.exit()
		else:
			args.genome=os.path.realpath(args.genome)
	else:
		args.genome="/mnt/genetica3/marius/pipeline_practicas_marius/hg19bundle"






	# checking if ped file exists

	if args.pedigree!=None:
		if not os.path.isfile(args.pedigree): 
			sys.stderr.write("ERROR: Panel file '%s' does not exist\n" %(args.pedigree))
			sys.exit()
		else:
			args.pedigree=os.path.realpath(args.pedigree)
			args.cvcf=True # combined genotyping
	else:
		args.pedigree="null"



	## basespace option not possible for bam files
	if args.skipMapping and args.basespace:
		sys.stderr.write("ERROR: Bam files can not be allocated in BaseSpace \n")
		sys.exit()





	# checking if correct csv method

	# if args.cnv-method!= None:
	# 	for i in args.cnv-method.split(","):
	# 		if i not in cnv_options:
	# 			sys.stderr.write("ERROR: At least one of the specified CNV methods is not correct: %s \n" %(args.cnv-method))
	# 			sys.exit()
	



	# checking and reading sample names within sample file (if provided)

	sample_names = list() #all project
	sample_namesT = list()	#tested samples
	cat=False

	print(args.samples)
	file_samples = list()
	if args.samples != "all":
		if os.path.isfile(args.samples):
			samplesFile = open(args.samples, "r")
			[file_samples.append((line.split(",")[1]).strip()) for line in samplesFile]
			file_samples = set(file_samples)
			if len(file_samples)==0:
				sys.stderr.write("ERROR: no samples taken from %s\n" %(args.samples))
				sys.exit()
			else:
				sys.stdout.write("\nUser-specified test samples: %s\n" %(", ".join(file_samples)))
		else:
			sys.stderr.write("ERROR: '%s' does not exist\n" %(args.samples))
			sys.exit()
	



	# when local samples, checking if folder exists and in case of fastq input check if need of concatenation.


	if not args.basespace:
		if not os.path.isdir(args.input): 
			sys.stderr.write("ERROR: Input folder '%s' does not exist\n" %(args.input))
			sys.exit()

		if not args.skipMapping:
			args.input=os.path.realpath(args.input)
			dnaids = [(os.path.basename(i)).split("_")[0].replace(".fastq.gz", "")  for i in glob(args.input+'/*.fastq.gz')]
			if len(dnaids)==0:
				sys.stderr.write("ERROR: Fastq files not found in '%s'\n" %(args.input))
			for dnaid in dnaids:
				sample_names.append(dnaid)			
				if args.samples == "all" or dnaid in file_samples:
					sample_namesT.append(dnaid)			
			# flag for fastq concatenation. 
			if args.analysis in ["cnv","all"]:
				fqs = Counter(sample_names).values()
			else: 
				fqs = Counter(sample_namesT).values()
			sample_namesT = set(sample_namesT)
			sample_names = set(sample_names)
			cat = not all(c == 2 for c in fqs)   # I could specify R1 and R2
		else:
			dnaids = [(os.path.basename(i)).split("_")[0].replace(".bam", "")  for i in glob(args.input+'/*.bam')]
			bai = [(os.path.basename(i)).split("_")[0].replace(".bai", "")  for i in glob(args.input+'/*.bai')]
			if len(dnaids) != len(set(dnaids)):
				sys.stderr.write("ERROR: more than one bam file per ID in folder '%s'\n" %(args.input))
				sys.exit()
			elif len(dnaids) != len(bai) or len(dnaids)==0:
				sys.stderr.write("ERROR: check input alignment directory '%s'. Different number of bam and bai files or not existing\n" %(args.input))
				sys.exit()
			for dnaid in dnaids:
				sample_names.append(dnaid)			
				if args.samples == "all" or dnaid in file_samples:
					sample_namesT.append(dnaid)	




	# define local dir for concatenated or basemount fastq files and run concatenation

	if (args.basespace or cat) and not args.skipMapping:
		inputDir = args.output + '/tmp_joinedFastq/'
		if not os.path.exists(inputDir):
	 		os.makedirs(inputDir)
	 	# run script to concatenate. 

 		sys.stdout.write("\n...............................................\n")
		sys.stdout.write("  CONCATENATION  (AND DOWNLOAD) OF FASTQ FILES     \n")
		sys.stdout.write("...............................................\n\n")

		myargs_desc = ["SCRIPT", "BASESPACE", "OUTPUT DIR", "PROJECT/INPUT FOLDER", "SAMPLE FILE", "ANALYSIS", "OUTPUT"]
		myargs = ["python "+utilitiesPath+"concatenation_1.py", str(args.basespace), inputDir, args.input, args.samples, args.analysis, args.output]
		[sys.stdout.write("%s: %s\n" %(myargs_desc[i], myargs[i])) for i in range(1,len(myargs)-1)]
		
		if args.local:	# LOCAL
		 	joinFastq(str(args.basespace), inputDir, args.input, args.samples, args.analysis, args.output)
		else: 
			job_name =  "FASTQret_"+run
			sys.stdout.write("JOB NAME: %s\n" %(job_name))
			jobid_conc = sbatch(job_name, args.output, ' '.join(myargs), time=4, mem=5, threads=args.threads, mail=args.mail, dep='')
			sys.stdout.write("JOB ID: %s\n" %(jobid_conc))
	else:
		inputDir = args.input



	if args.basespace:
		dnaids = set([(os.path.basename(i)).split("_")[0].replace(".fastq.gz", "")  for i in glob(inputDir+'/*.fastq.gz')])
		if len(dnaids)==0:
			sys.stderr.write("ERROR: Fastq files not found in '%s'\n" %(args.input))
		for dnaid in dnaids:
			sample_names.append(dnaid)			
			if args.samples == "all" or dnaid in file_samples:
				sample_namesT.append(dnaid)	



	if args.samples != "all":
		if len(list(set(file_samples).intersection(sample_namesT))) != len(file_samples):
			sys.stdout.write("\nWARNING: This samples could not be found inside the input directory: \n%s\n\n" %(",".join(set(file_samples) - set(sample_namesT))))
		if len(set(sample_namesT))==0:
			sys.stderr.write("ERROR: no fastq.gz/bam files to analyse. \n")
			sys.exit()
		else:
			sys.stdout.write("\nSamples that will be analysed: %s\n" %(",".join(set(sample_namesT))))
	else:
		sys.stdout.write("\nSamples that will be analysed: %s\n" %(",".join(set(sample_namesT))))


	sample_namesT = set(sample_namesT)
	sample_names = set(sample_names)


	# Analysing individual examples if mapping or SNP are specified


	jobid_list=[]
	jobid_list_snp=[]
	depJobs=''
	
	if args.skipMapping == False or args.analysis in ["snp","all"]: # mapping or snps

		sys.stdout.write("\n......................................................\n")
		sys.stdout.write("  RUNNING MAPPING OR/AND SNV CALLING FOR INDIVIDUAL SAMPLES\n")
		sys.stdout.write("......................................................\n")


		sample_string = ",".join(sample_names)
		sys.stdout.write("Samples: %s\n" %(sample_string))


 		for sample_name in sample_names: # four options: (1) no sample analysis; (2) just mapping; (3) just snp; (4) mapping + snp

			if args.analysis in ["cnv","all"] or sample_name in sample_namesT: # just analyse if contained in samples_files or we are running CNVs	

				if args.skipMapping and sample_name in sample_namesT:
					sampleAnalysis = "snp" 

				elif args.analysis in ["snp","all"] and sample_name in sample_namesT:
					sampleAnalysis = "mapping_snp" 
					
				else:
					sampleAnalysis = "mapping"

 				### MAPPING AND GENOTYPING
				sys.stdout.write("\nAnalysing individual sample '%s' (%s) with arguments:\n" %(sample_name, sampleAnalysis))
				myargs_desc = ["SCRIPT", "INPUT FOLDER", "OUTPUT DIR", "SAMPLE LABEL", "N THREADS", "RUN LABEL", "PANEL BED FILE", "ANALYSIS TYPE", "COMBINED VCF", "SKIPPING MAPPING STEP", "GENOME BUNDLE", "LOCAL", "PATHOLOGY", "INTERVALS"]
				myargs = [utilitiesPath+"SNV_pipeline.sh", inputDir, args.output, sample_name, str(args.threads), run, args.panel, sampleAnalysis, str(args.cvcf), str(args.skipMapping), args.genome, str(args.local), str(args.pathology), str(args.intervals)]
				[sys.stdout.write("%s: %s\n" %(myargs_desc[i], myargs[i])) for i in range(1,len(myargs))]
				job_name =  sampleAnalysis+"_"+sample_name
				sys.stdout.write("JOB NAME: %s\n" %(job_name))

				start_time = time.time()

				if args.local:	# LOCAL
					stdout_f = open(args.output+"/"+job_name+".out", 'w')
					stderr_f = open(args.output+"/"+job_name+".err", 'w')
				 	subprocess.call(myargs, stdout= stdout_f, stderr = stderr_f)
				
				elif sampleAnalysis=="mapping": # SBATCH AND SAVE JOB IDS FOR MAPPING SAMPLES
					jobid=sbatch(job_name, args.output, ' '.join(myargs), time=400, mem=5, threads=args.threads, mail=args.mail, dep='')
					print(jobid)
					jobid_list.append(jobid)
					sys.stdout.write("JOB ID: %s\n" %(jobid))

				else: # SBATCH AND SAVE JOB IDS FOR MAPPING AND SNP CALLING SAMPLES
					jobid_snp=sbatch(job_name, args.output, ' '.join(myargs), time=400, mem=5, threads=args.threads, mail=args.mail, dep='')
					jobid_list_snp.append(jobid_snp)
					sys.stdout.write("JOB ID: %s\n" %(jobid_snp))

				print("TOTAL SECONDS: %s" % (time.time() - start_time))




	# Joint Genotyping for multiple sample analysis or families  

	if args.cvcf:

		sys.stdout.write("\n...............................\n")
		sys.stdout.write("  RUNNING COMBINED GENOTYPING   \n")
		sys.stdout.write("...............................\n\n")


		sys.stdout.write("Combining genotyping with arguments:\n")
		myargs_desc = ["SCRIPT", "OUTPUT DIR", "RUN LABEL", "NUMBER THREADS", "GENOME BUNDLE","PEDIGREE","LOCAL","PATHOLOGY"]
		myargs_cvcf = [utilitiesPath+"combinedGenotyping.sh", args.output, run, str(args.threads), args.genome, args.pedigree, str(args.local), str(args.pathology)]

		[sys.stdout.write("%s: %s\n" %(myargs_desc[i], myargs_cvcf[i])) for i in range(1,len(myargs_cvcf))]
		sys.stdout.write("FILE WITH SAMPLE NAMES: %s/haplotype_caller_gvcf_data/my_list_of_gvcfs_files_to_combine_%s.list\n" %(args.output, run))

		job_name = "cVCF_"+run
		sys.stdout.write("JOB NAME: %s\n" %(job_name))


		start_time = time.time()

		if args.local:
			stdout_f = open(args.output+"/"+job_name+".out", 'w')
			stderr_f = open(args.output+"/"+job_name+".err", 'w')
			subprocess.call(myargs_cvcf, stdout= stdout_f, stderr = stderr_f)
		
		else:
			depJobs = ':'.join(jobid_list_snp)  # cVCF just dependes on samples jobs id that have been mapped and snp called.
			jobid2=sbatch(job_name, args.output, ' '.join(myargs_cvcf), time=400, mem=10, threads=args.threads, mail=args.mail, dep=depJobs)
			sys.stdout.write("JOB ID: %s\n" %(jobid2))

		print("TOTAL SECONDS: %s" % (time.time() - start_time))






	# Copy Number Variant analysis

	if args.analysis=="cnv" or args.analysis=="all":


		sys.stdout.write("\n.............................................\n")
		sys.stdout.write("  RUNNING CNV DETECTION FOR PROVIDED SAMPLES \n")
		sys.stdout.write(".............................................\n\n")

		if args.skipMapping:
			bamF = inputDir
		else:
			bamF = args.output+"/applied_bqsr_data"
		sys.stdout.write("Copy number variants calling with arguments:\n")
		myargs_desc = ["SCRIPT", "OUTPUT DIR", "BAM DIRECTORY", "SAMPLES FILE", "RUN LABEL", "NUMBER THREADS", "PANEL BED FILE", "WINDOW"]
		myargs_cnv = [utilitiesPath+"CNVdetection.sh", args.output, bamF,  args.samples, run, str(args.threads), args.panel, "no", utilitiesPath, str(args.local)]

		[sys.stdout.write("%s: %s\n" %(myargs_desc[i], myargs_cnv[i])) for i in range(1,len(myargs_cnv)-2)]
		sys.stdout.write("BAM FILES:\n%s" %("\n".join(glob(bamF + '/*.bam'))))

		job_name = "CNV_"+run
		sys.stdout.write("\nJOB NAME: %s\n" %(job_name))


		start_time = time.time()

		if args.local:	
			stdout_f = open(args.output+"/"+job_name+".out", 'w')
			stderr_f = open(args.output+"/"+job_name+".err", 'w')
			subprocess.call(myargs_cnv, stdout= stdout_f, stderr = stderr_f)
		
		else:
			if len(jobid_list+jobid_list_snp)==0: # in case all samples are mapped, cnv calling don't have dependent jobs.
				depJobs==''
			else:
				depJobs = ':'.join(jobid_list+jobid_list_snp)
			 # cnv depends on all samples "mapped" and "mapped + snp"
			jobid3=sbatch(job_name, args.output, ' '.join(myargs_cnv), time=400, mem=10, threads=args.threads, mail=args.mail, dep=depJobs)
			sys.stdout.write("JOB ID: %s\n" %(jobid3))
			sys.stdout.write("DEPENDENT JOBS: %s\n" %(depJobs))
	
		print("TOTAL SECONDS: %s" % (time.time() - start_time))


	sys.stdout.close()
	sys.stderr.close()




	# mycmdcopySW = 'cat %s >> %s' %(args.output+"/software"+run+".txt", stdoutAll)
	# subprocess.Popen(mycmdcopySW, shell=True)



if __name__ == "__main__":
	main()

