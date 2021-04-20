#!/usr/bin/env python

#######################################################################
# Script: Run FJD analysis for samples stored in basespace or locally
# Author: Lorena de la Fuente 
# Date: 24-01-2019
#######################################################################

import sys
import os
import re
from glob import glob
import subprocess
import argparse
from collections import Counter
import datetime
import time
from collections import defaultdict 

now = (datetime.datetime.now()).strftime('%Y_%m_%d_%H_%M_%S')
print("DATE:"+now)
softwarePath = os.path.dirname(os.path.realpath(__file__))+"/"
pipelinesPath =  os.path.dirname(os.path.realpath(__file__))+"/pipelines/" 
tasksPath =  os.path.dirname(os.path.realpath(__file__))+"/tasks/" 



def sbatch(job_name, folder_out, command, mem=4, time=400, threads=1, mail=None, dep='', wait = '', queue=None, typedep=None):

	if dep != '':
		if typedep!=None:
			dep = '--dependency=afterany:{} --kill-on-invalid-dep=yes '.format(dep)
		else:
			dep = '--dependency=afterok:{} --kill-on-invalid-dep=yes '.format(dep)

	if mail!=None:
		mailc = "--mail-user={} --mail-type=FAIL".format(mail)
	else:
		mailc = ''

	if wait==True:
		wait = '--wait'

	if queue==True:
 		partition='--partition=fastbioinfo'
 	else:
 		partition='--partition=bioinfo'


	sbatch_command = "sbatch -J {} -o {}/{}.out -e {}/{}.err {} -t {}:00:00 --account=bioinfo_serv {} --mem-per-cpu={}gb --cpus-per-task={} {} {} {}".format(job_name, folder_out, job_name, folder_out, job_name, mailc, time, partition, mem, threads, dep, wait, command)
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
	parser.add_argument('-a', '--analysis', help='\t\tType of analysis to run', required=False, choices={"mapping", "snv", "cnv", "all"}, default="snv")
	parser.add_argument('-c', '--cvcf', help='\t\tCombined genotyping. Number of samples must be higher than 2', required=False, action='store_true')
	parser.add_argument('-p', '--panel', help='\t\tBed with panel regions. Mandatory if CNV analysis', required=False)
	parser.add_argument('-b', '--basespace', help='\t\tTake samples from Basespace', required=False, action='store_true')
	parser.add_argument('-l', '--local', help='\t\tRun in local', required=False, action='store_true')
	parser.add_argument('-t', '--threads', help='\t\tNumber of threads for mapping', type=int, required=False, default=1)
	parser.add_argument('-M', '--memory', help='\t\tNumber of GBs for Java GATK and VEP running - Total GBs', type=int, required=False, default=20)
	parser.add_argument('-f', '--genefilter', help='\t\tGene list to filter SNVs', required=False, default=False)
	parser.add_argument('-T', '--time', help='\t\tNumber of hours', type=int, required=False, default=1000)
	parser.add_argument('-m', '--mail', help='\t\tMail account', required=False)
	parser.add_argument('-S', '--single' , help='\t\tConserve sample VCFs when running combined genotyping (CVCF)', required=False, action='store_true')
	parser.add_argument('-k', '--skipMapping', help='\t\tOption to skip mapping. Input folder must be fastq folder and Output folder the global output folder', action='store_true')
	parser.add_argument('-g', '--genome', help='\t\tLocal directory for genome reference file', required=False)
	parser.add_argument('-e', '--pedigree', help='\t\tPedigree file. Combined genotyping will be run automatically', required=False)
	parser.add_argument('-v', '--version', help="Display program version number.", action='version', version='1.0')
	parser.add_argument('-P', '--pathology', help="Disease group name: resp, digest, skin, musc, gu, pregnancy, perinatal, cong_mal, clinic, infectious, other, neoplasm, blood, endoc, mental, CNS, eye, ear, circ,healthy.", default='healthy')
	parser.add_argument('-I', '--intervals', help="Specified if genomic intervals (panel) over which to operate during SNV calling are used.", action='store_true')
	parser.add_argument('-C', '--cnv-method', help="Method to call copy number variants (CNVs). Comma-separated list of methods. Methods: ED, CN, C2, PM.", required=False, default="ED,CN,C2,PM")
	parser.add_argument('-w', '--window', help="125 window size for CNV intervals", required=False, action='store_true')
	parser.add_argument('-B', '--remove_bam', help="Remove bam files", required=False, action='store_true')
	parser.add_argument('-G', '--padding', help='\t\tPadding length', type=int, required=False, default=1000)
	parser.add_argument('-d', '--duplicates', help="Data is not dedupped. By default, duplicates are marked", required=False, action='store_true')
	parser.add_argument('-u', '--basemountuser', help="Specify alternative basemount user", required=False, default=False)
	parser.add_argument('-Q', '--qcnvthreshold', help='\t\tMinimum read coverage for 90percent of targets for CNV analysis.', required=False, default=40)
	parser.add_argument('-q', '--queue', help='\t\tUse the fastqueue.', required=False, action="store_true")





	fjd_start_time = time.time()
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




	## basespace option not possible for bam files
	
	if args.skipMapping and args.basespace:
		sys.stderr.write("ERROR: Bam files can not be retrieved from BaseSpace \n")
		sys.exit()



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
			sys.exit()

		if args.analysis in ["cnv","all"]:
			sys.stderr.write("ERROR: Bed file of target regions is mandatory for CNV analysis.\n")
			sys.exit()




	# checking if specified bundle exists

	if args.genome!=None:
		if not os.path.isdir(args.genome): 
			sys.stderr.write("ERROR: Genome bundle folder '%s' does not exist\n" %(args.genome))
			sys.exit()
		else:
			args.genome=os.path.realpath(args.genome)
	else:
		if args.local:
			args.genome="/mnt/genetica3/marius/pipeline_practicas_marius/hg19bundle/ucsc.hg19.fasta"
		else:
			args.genome="/home/proyectos/bioinfo/references/hg19/ucsc.hg19.fasta"





	# checking if ped file exists

	if args.pedigree!=None:
		if not os.path.isfile(args.pedigree): 
			sys.stderr.write("ERROR: Panel file '%s' does not exist\n" %(args.pedigree))
			sys.exit()
	else:
		args.pedigree="null"



	# checking if gene list exits
	if args.genefilter!=False:
		if not os.path.isfile(args.genefilter): 
			sys.stderr.write("ERROR: File listing genes of interest ('%s') does not exist\n" %(args.genefilter))
			sys.exit()



	# checking if bam folder is empty

	removebam=False

	if args.skipMapping:
		bamF = args.input
	else:
		bamF = args.output+"/bams"
	
	if not args.skipMapping and args.remove_bam:
		if args.analysis != ["cnv","all"]:
			removebam="single"
		else:
			if not os.path.isdir(bamF) or not os.listdir(bamF):
				removebam="folder"





	# checking if correct csv method
	
	cnv_options=["ED","CN","C2","PM"]
	if args.cnv_method!="ED,CN,C2,PM":
		for i in args.cnv_method.split(","):
			if i.strip() not in cnv_options:
				sys.stderr.write("ERROR: At least one of the specified CNV methods is not correct: %s \n" %(args.cnv_method))
				sys.exit()
	



	# checking and reading sample names within sample file (if provided)

	sample_names = list() #all project
	sample_namesT = list()	#tested samples
	cat=False

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
				sys.stdout.write("\n..................................\n")
				sys.stdout.write("User-specified test samples:\n%s" %(", ".join(file_samples)))
				sys.stdout.write("\n..................................\n")

		else:
			sys.stderr.write("ERROR: '%s' does not exist\n" %(args.samples))
			sys.exit()
	



	# when local samples, checking if folder exists and in case of fastq input check if need of concatenation.


	if not args.basespace:

		if not os.path.isdir(args.input): 
			sys.stderr.write("ERROR: Input folder '%s' does not exist\n" %(args.input))
			sys.exit()
		else:
			args.input=os.path.realpath(args.input)

		if not args.skipMapping:
			# if fasta provided, the name would be the string before "_"
			files = glob(args.input+'/*.f*q.gz')

			dnaids=[]
			dnaidsR1=defaultdict(int)
			dnaidsR2=defaultdict(int)

			dnaids = set(sorted([(os.path.basename(i)).split("_")[0].replace(".fastq.gz", "").replace(".fq.gz", "")  for i in files]))
			for w in [(os.path.basename(i)).split("_")[0].replace(".fastq.gz", "").replace(".fq.gz", "") for i in files if '_R1' in i]: dnaidsR1[w] +=1
			for w in [(os.path.basename(i)).split("_")[0].replace(".fastq.gz", "").replace(".fq.gz", "") for i in files if '_R2' in i]: dnaidsR2[w] +=1

			if len(dnaids)==0:
				sys.stderr.write("ERROR: Fastq files not found in '%s'\n" %(args.input))
				sys.exit()
			for dnaid in dnaids:
				sample_names.append(dnaid)			
				if args.samples == "all" or dnaid in file_samples:
					sample_namesT.append(dnaid)
				if args.analysis in ["cnv","all"] or dnaid in file_samples or args.samples=="all":
					if dnaid in dnaidsR1.keys() and dnaid in dnaidsR2.keys():
						if dnaidsR1[dnaid]!= dnaidsR2[dnaid]:
							sys.stderr.write("ERROR: check fastq input directory '%s'. Different number of files for R1 and R2\n" %(args.input))
							sys.exit()
						elif dnaidsR1[dnaid]>1 or dnaidsR2[dnaid]>1:
							cat=True
					else:
						sys.stderr.write("ERROR: check fastq input directory '%s'. No files for R1 or R2.\n" %(args.input))
						sys.exit()
		else:
			dnaids = [(os.path.basename(i)).replace(".bam", "")  for i in glob(args.input+'/*.bam')]
			bai = [(os.path.basename(i)).replace(".bai", "")  for i in glob(args.input+'/*.bai')]
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
	
	else:

		# Check project name and if exists list all samples associated to that project

		sys.stdout.write("\nChecking if project exists...\n")

		projectfile = args.output+'/'+run+'_projects.txt'
		datasetfile = args.output+'/'+args.input+'_'+run+'_'+'datasets.txt'

		projName = "--project-name="+args.input
		stdoutFile = "--stdout="+datasetfile

		if args.basemountuser==False:
			
			config = ""

		else:

			config = "--config "+args.basemountuser


		# bs list projects

		cmd="/home/proyectos/bioinfo/software/bs list project "+config+ " > "+ projectfile
		subprocess.call(cmd, shell=True)


		with open(projectfile) as myfile:
			if args.input not in myfile.read().split():
				sys.stderr.write("ERROR: Project '%s' does not exist in basespace\n" %(args.input))
			else:
				os.remove(projectfile)


		### bs list dnaids (samples)

		cmdSamples = " ".join(["/home/proyectos/bioinfo/software/bs",  "list", "biosample", "--sort-by=BioSampleName", config, projName,  stdoutFile])
		subprocess.call(cmdSamples, shell=True)


		# if project samples, we save the entire biosample name but check with the input sample file by looking at the first 7 characters 

		dnaids_dicc = {}
		with open(datasetfile) as myfile:
			lines = myfile.read().splitlines()
			lines = lines[3:len(lines)-1]
			for line in lines:
				biosampleName = line.split()[1]
				dnaid = biosampleName[0:7]
				dnaids_dicc[dnaid] = biosampleName
		for dnaid,biosampleName in dnaids_dicc.iteritems():
			sample_names.append(biosampleName)			
			if args.samples == "all" or dnaid in file_samples:
				sample_namesT.append(biosampleName)
		dnaids = dnaids_dicc.keys()	


	# Sample summary

	if args.samples != "all":
		if len(list(set(file_samples).intersection(sample_namesT))) != len(file_samples):
			sys.stdout.write("\nWARNING - This samples could not be found inside the input directory: \n%s\n\n" %(",".join(set(file_samples) - set(dnaids))))
		if len(set(sample_namesT))==0:
			sys.stderr.write("ERROR: no fastq.gz/bam files to analyse. \n")
			sys.exit()


	sample_namesT = set(sample_namesT)
	sample_names = set(sample_names)

	sys.stdout.write("\n..................................\n")
	sys.stdout.write("Samples that will be analysed: \n%s" %(", ".join(set(sample_namesT))))
	sys.stdout.write("\n..................................\n")


	if len(sample_names - sample_namesT) != 0 and args.analysis in ["cnv","all"] :
		sys.stdout.write("\n..................................\n")
		sys.stdout.write("Control Samples: \n%s" %(",".join(set(sample_names))))
		sys.stdout.write("\n..................................\n")


	if args.basespace or cat:
		if args.local:
			inputDir = args.output + '/tmp_joinedFastq/'
		else:	
			inputDir = '/scratch/' + os.environ["USER"] + '/' + run + '/'
		if not os.path.exists(inputDir):
			os.makedirs(inputDir)
	else: 
		inputDir = args.input










	#*********** 1. Analysing individual examples (mapping or/and SNV calling)

	jobid_list_mapProc=[]
	jobid_list_snpAnnot=[]
	depJobs=''
	

	if args.skipMapping == False or args.analysis in ["snv","all"]: # mapping or snps

		sys.stdout.write("\n........................................................\n")
		sys.stdout.write("RUNNING MAPPING OR/AND SNV CALLING FOR INDIVIDUAL SAMPLES\n")
		sys.stdout.write("........................................................\n")


		for sample_name in sample_names: # four options: (1) no sample analysis; (2) just mapping; (3) just snp; (4) mapping + snp

			sampleAnalysis = []


			if args.analysis in ["cnv","all"] or sample_name in sample_namesT: # sample is analyzed if contained in samples_files or we are running CNVs (need of controls). Depending on that we assign the type of analysis to the sample

				if args.skipMapping:
					
					if sample_name in sample_namesT:
						sampleAnalysis.append("snv") 

				elif args.analysis in ["snv","all"] and sample_name in sample_namesT:
					sampleAnalysis=["mapping","snv"]
	
				else:
					sampleAnalysis = ["mapping"]


					
				if len(sampleAnalysis)!=0:

					sys.stdout.write("\nAnalysing individual sample '%s' (%s) with arguments:\n\n" %(sample_name, ":".join(sampleAnalysis)))
					
					# print all arguments in stdout
					myargs_desc = ["INPUT FOLDER", "OUTPUT DIR", "SAMPLE LABEL", "N THREADS", "RUN LABEL", "PANEL BED FILE", "BASESPACE DOWNLOAD", "CONCATENATION", "INPUT DIRECTORY", "ANALYSIS TYPE", "COMBINED VCF", "SKIPPING MAPPING STEP", "GENOME BUNDLE", "LOCAL", "PATHOLOGY", "INTERVALS", "KEEP DUPLICATES", "REMOVE BAM FILES", "GENE LIST", "PADDING", "BASESPACE USER", "SOFTPATH"]
					myargs = [args.input, args.output, sample_name, str(args.threads), run, args.panel, str(args.basespace), str(cat), inputDir, sampleAnalysis, str(args.cvcf), str(args.skipMapping), args.genome, str(args.local), str(args.pathology), str(args.intervals), str(args.duplicates), str(removebam), str(args.genefilter), str(args.padding), str(args.basemountuser) ,softwarePath]
					[sys.stdout.write("%s: %s\n" %(myargs_desc[i], myargs[i])) for i in range(1,len(myargs)-1)]
					

					myargs_pipe1 = [pipelinesPath+"pipeline1_downloadMapping.sh", args.input, args.output, sample_name, str(args.threads), run, str(args.basespace), str(cat), inputDir,  args.genome, str(args.local), str(args.basemountuser), softwarePath]
					myargs_pipe1_2 = [tasksPath+"BAMpreprocessing.sh",  str(args.local), run, args.output, sample_name, str(args.duplicates), args.genome, str(args.memory)]
					myargs_pipe2 = [pipelinesPath+"pipeline2_QCbamSNVCallingFiltering.sh", args.input, args.output, sample_name, str(args.memory), run, args.panel, str(args.cvcf), str(args.skipMapping), args.genome, str(args.local), str(args.intervals), str(removebam),  str(args.padding), softwarePath]
					myargs_pipe4 = [pipelinesPath+"pipeline4_LohAnnotationOutput.sh", args.output, sample_name, str(4), run, args.panel, str(args.cvcf), args.genome, str(args.local), args.pathology, str(args.genefilter), softwarePath]


					start_time = time.time()
					procid = ""


					if "mapping" in sampleAnalysis: # SBATCH AND SAVE JOB IDS FOR MAPPING SAMPLES
						
						jobname_pipe1 =  "mapping"+"_"+sample_name
						mapid=sbatch(jobname_pipe1, args.output, ' '.join(myargs_pipe1), time=args.time, threads=args.threads, mail=args.mail, dep='') # mmapping alone: just to take advantage of threads
						jobid_list_mapProc.append(mapid)
						sys.stdout.write("JOB ID MAPPING: %s\n" %(mapid))

						jobname_pipe2 =  "preprocessing_"+sample_name
						procid=sbatch(jobname_pipe2, args.output, ' '.join(myargs_pipe1_2), time=args.time, threads=args.memory/4, mail=args.mail, dep=mapid) # bam processing + snp calling + snp filtering
						jobid_list_mapProc.append(procid)
						sys.stdout.write("JOB ID PRE-PROCESSING: %s\n" %(procid))


					if "snv" in sampleAnalysis: # SBATCH AND SAVE JOB IDS FOR SNP CALLING SAMPLES


						jobname_pipe2 =  "snv_"+sample_name
						snpid=sbatch(jobname_pipe2, args.output, ' '.join(myargs_pipe2), time=args.time, threads=args.memory/4, mail=args.mail, dep=procid) # bam processing + snp calling + snp filtering
						jobid_list_snpAnnot.append(snpid)
						sys.stdout.write("JOB ID SNV CALLING: %s\n" %(snpid))


						if not args.cvcf:

							jobname_pipe4 =  "annotation_"+sample_name
							annotid=sbatch(jobname_pipe4, args.output, ' '.join(myargs_pipe4), time=args.time, threads=args.memory/4, mail=args.mail, dep=snpid) # vcf annotation and complementary analysis to generate full output
							jobid_list_snpAnnot.append(annotid)
							sys.stdout.write("JOB ID SNV ANNOTATION: %s\n" %(annotid))



					print("TOTAL SECONDS: %s\n" % (time.time() - start_time))






	#*********** 2. Joint Genotyping for multiple sample analysis or families  


	cvcfjobs_list = []

	if args.cvcf:

		sys.stdout.write("\n...............................\n")
		sys.stdout.write("  RUNNING COMBINED GENOTYPING   \n")
		sys.stdout.write("...............................\n\n")

		myargs_pipe3 = [pipelinesPath+"pipeline3_combinedGenotyping.sh", args.output, str(args.local), run, args.genome, str(args.cvcf), args.pedigree, str(args.memory), softwarePath]
		jobname_pipe3 =  "combinedGenotyping"+"_"+run
		depJobs = ':'.join(jobid_list_snpAnnot)  	
		combid=sbatch(jobname_pipe3, args.output, ' '.join(myargs_pipe3), time=args.time, mem=4, threads=1, mail=args.mail, dep=depJobs) # mmapping alone: just to take advantage of threads
		cvcfjobs_list.append(combid)
		sys.stdout.write("JOB ID COMBINED GENOTYPING: %s\n" %(combid))



		# Annotation and processing of output 
		
		myargs_pipe4[3] = run
		job_name="annotation_"+run
		jobid=sbatch(job_name, args.output, ' '.join(myargs_pipe4), time=args.time, mem=args.memory/4, threads=4, mail=args.mail, dep=combid) 
		cvcfjobs_list.append(jobid)
		sys.stdout.write("JOB ID SNV ANNOTATION: %s\n" %(jobid))

		



	#*********** 3. If everything went well we should copy results to MAF folder!!!!!!!!!!!




	#*********** 4. Copy number variation for panel and exomes

	vdepCVjobs_list = []
	jobidQC_list = []


	if args.analysis=="cnv" or args.analysis=="all":

		start_time = time.time()

		sys.stdout.write("\n.............................................\n")
		sys.stdout.write("RUNNING CNV DETECTION FOR PROVIDED SAMPLES \n")
		sys.stdout.write(".............................................\n\n")

		method = "QC,"+args.cnv_method+",MA"

		# write arguments for users 
		sys.stdout.write("Copy number variants calling with arguments:\n")
		myargs_desc = [ "OUTPUT DIR", "BAM DIRECTORY", "SAMPLES FILE", "RUN LABEL", "NUMBER THREADS", "PANEL BED FILE", "WINDOW", "TASK FOLDER", "LOCAL", "METHODS", "GENE LIST", "GENOME"]
		myargs_cnv = [ args.output, bamF,  args.samples, run, str(args.threads), args.panel, str(args.window), tasksPath, str(args.local), method, str(args.genefilter), str(args.genome)]
		[sys.stdout.write("%s: %s\n" %(myargs_desc[i], myargs_cnv[i])) for i in range(1,len(myargs_cnv))]
		sys.stdout.write("BAM FILES:\n%s" %("\n".join(glob(bamF + '/*.bam'))))


		# run CNV-associated tasks after dependent jobs

		depJobs_list = jobid_list_mapProc

		for method in method.split(","):
			myargs_cnv = [pipelinesPath+"pipeline10_CNVcalling.sh", args.input, args.output, str(args.skipMapping), args.samples, run, str(args.threads), args.panel, str(args.window), tasksPath, str(args.local), method,   str(args.genefilter), str(args.genome), str(args.qcnvthreshold)]
			job_name = method+"_CNV_"+run

			if method == "QC":
				sys.stdout.write("\n\n#########  JOB FOR QUALITY CONTROL OF PANEL FILE AND SAMPLES \n")
				depJobs = ':'.join(depJobs_list)
				jobidQC=sbatch(job_name, args.output, ' '.join(myargs_cnv), time=args.time, mem=args.memory, threads=2, mail=args.mail, dep=depJobs, typedep="any")
				jobidQC_list.append(jobidQC)
				sys.stdout.write("JOB ID: %s\n" %(jobidQC))
				sys.stdout.write("DEPENDENT JOBS: %s\n" %(depJobs))
			elif method == "MA":
				sys.stdout.write("\n#########  COMBINING CNV RESULTS FROM ALTERNATIVE METHODS \n")
				depJobs = ':'.join(vdepCVjobs_list+[jobidQC])
				jobidMA=sbatch(job_name, args.output, ' '.join(myargs_cnv), time=args.time, mem=args.memory, threads=2, mail=args.mail, dep=depJobs, typedep="any")
				vdepCVjobs_list.append(jobidMA)
				sys.stdout.write("JOB ID: %s\n" %(jobidMA))
				sys.stdout.write("DEPENDENT JOBS: %s\n" %(depJobs))
			else:
				sys.stdout.write("\n#########  CNV CALLING FOR %s\n" %(method))
				jobidME=sbatch(job_name, args.output, ' '.join(myargs_cnv), time=args.time, mem=args.memory, threads=2, mail=args.mail, dep=jobidQC)
				vdepCVjobs_list.append(jobidME)
				sys.stdout.write("JOB ID: %s\n" %(jobidME))
				sys.stdout.write("DEPENDENT JOBS: %s\n" %(jobidQC))

			print(job_name)

		print("TOTAL SECONDS: %s" % (time.time() - start_time))





	#*********** 5. Remove empty directories


	myargs_remove = [tasksPath+"removeDirs.sh", args.output, str(removebam), bamF]


	sys.stdout.write("\n................................\n")
	sys.stdout.write("REMOVING TMP FOLDERS AND FILES \n")
	sys.stdout.write("................................\n\n")


	depJobs_list = jobid_list_mapProc + jobid_list_snpAnnot + jobidQC_list + vdepCVjobs_list + cvcfjobs_list
	depJobs = ':'.join(depJobs_list)
	job_name = "removeDirs_"+run
	jobidREMOVE=sbatch(job_name, args.output, ' '.join(myargs_remove), time=args.time, mem=1, threads=1, mail=args.mail, dep=depJobs, typedep="any")
	depJobs_list.append(jobidREMOVE)
	sys.stdout.write("JOB ID: %s\n" %(jobidREMOVE))
	sys.stdout.write("DEPENDENT JOBS: %s\n" %(depJobs))







	#*********** 6. Running summary


	# create a new job which analysed what happend with all the jobs ran during these task

	sys.stdout.write("\n...........................\n")
	sys.stdout.write("SLURM JOB STATUS SUMMARY \n")
	sys.stdout.write("...........................\n\n")

	depJobs = ':'.join(depJobs_list)
	job_name = "slurmSummary_"+run

	myargs_summary = [tasksPath+"checkproject.sh", args.output, depJobs, run]

	jobidSUM=sbatch(job_name, args.output, ' '.join(myargs_summary), time=args.time, mem=1, threads=1, mail=args.mail, dep=depJobs, typedep="any")
	sys.stdout.write("JOB ID: %s\n" %(jobidSUM))
	sys.stdout.write("DEPENDENT JOBS: %s\n" %(depJobs))










	sys.stdout.close()
	sys.stderr.close()


if __name__ == "__main__":
	main()

