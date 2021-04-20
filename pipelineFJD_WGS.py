#!/usr/bin/env python

######################################################
# Script: Run FJD analysis for WGS local samples
# Author: Lorena de la Fuente 
# Date: 16-03-2021
######################################################


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
utilitiesPath =  os.path.dirname(os.path.realpath(__file__))+"/utilities/" 
pipelinesPath =  os.path.dirname(os.path.realpath(__file__))+"/pipelines/" 
tasksPath =  os.path.dirname(os.path.realpath(__file__))+"/tasks/" 


def sbatch(job_name, folder_out, command, mem=5, time=400, threads=1, mail=None, dep='', wait = '', typedep=None):

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


	sbatch_command = "sbatch -J {} -o {}/{}.out -e {}/{}.err {} -t {}:00:00 --account=bioinfo_serv --partition=bioinfo --mem-per-cpu={}gb --cpus-per-task={} {} {} {}".format(job_name, folder_out, job_name, folder_out, job_name, mailc, time, mem, threads, dep, wait, command)
	sbatch_response = subprocess.check_output(sbatch_command, shell=True)
	job_id = sbatch_response.split(' ')[-1].strip()
	return job_id
			




def main():

	#arguments
	parser = argparse.ArgumentParser(description="FJD WGS analysis")
	parser.add_argument('-i', '--input', help='\t\t Local directory with fastq.gz/bam input files', required=True)
	parser.add_argument('-o', '--output', dest="output",help='\t\tOutput directory name.', required=True)
	parser.add_argument('-s', '--samples' , help='\t\tCSV file listing sample names within second column', required=False, default="all")
	parser.add_argument('-n', '--name', help='\t\tName for job. RECOMMENDED', required=False)
	parser.add_argument('-a', '--analysis', help='\t\tType of analysis to run', required=False, choices={"mapping", "snv", "all"}, default="snv")
	parser.add_argument('-c', '--cvcf', help='\t\tCombined genotyping. Number of samples must be higher than 2', required=False, action='store_true')
	parser.add_argument('-l', '--local', help='\t\tRun in local', required=False, action='store_true')
	parser.add_argument('-t', '--threads', help='\t\tNumber of threads', type=int, required=False, default=4)
	parser.add_argument('-M', '--memory', help='\t\tNumber of GBs for VEP running', type=int, required=False, default=5)
	parser.add_argument('-T', '--time', help='\t\tNumber of hours', type=int, required=False, default=1000)
	parser.add_argument('-m', '--mail', help='\t\tMail account', required=False)
	parser.add_argument('-S', '--single' , help='\t\tConserve sample VCFs when running combined genotyping (CVCF)', required=False, action='store_true')
	parser.add_argument('-k', '--skipMapping', help='\t\tOption to skip mapping. Input folder must be fastq folder and Output folder the global output folder', action='store_true')
	parser.add_argument('-g', '--genome', help='\t\tLocal directory for genome reference file', required=False)
	parser.add_argument('-e', '--pedigree', help='\t\tPedigree file. Combined genotyping will be run automatically', required=False)
	parser.add_argument('-P', '--pathology', help="Disease group name: resp, digest, skin, musc, gu, pregnancy, perinatal, cong_mal, clinic, infectious, other, neoplasm, blood, endoc, mental, CNS, eye, ear, circ,healthy.", default='healthy')
	parser.add_argument('-I', '--intervals', help="Specified if genomic intervals (panel) over which to operate during SNV calling are used.", action='store_true')
	parser.add_argument('-D', '--depth', help="Skip filtering of samples during CNV calling by read coverage.", required=False, action='store_true')
	parser.add_argument('-B', '--remove_bam', help="Remove bam files", required=False, action='store_true')
	parser.add_argument('-f', '--genefilter', help='\t\tRestricted analysis to genes listed in file', required=False, default=False)
	parser.add_argument('-v', '--version', help="Display program version number.", action='version', version='2.0')


	fjd_start_time = time.time()
	args = parser.parse_args()


	# label for current run

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

	sys.stdout.write("\nFJD ANALYSIS FOR WHOLE GENOME SEQUENCING\n")
	sys.stdout.write("\nChecking arguments...\n")



	# checking output directory

	if not os.path.isdir(args.output): 
		sys.stderr.write("ERROR: Output folder '%s' does not exist\n" %(args.output))
		sys.exit()
	else:
		args.output=os.path.realpath(args.output)



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
			intervals_files = glob('/home/proyectos/bioinfo/references/b37_scattered_wgs_intervals/scattered_wgs_intervals/hg19_b37_scattered_wgs_intervals_scatter*_list.bed')



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




	# checking and reading sample names within sample file (if provided)

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

	if not os.path.isdir(args.input): 
		sys.stderr.write("ERROR: Input folder '%s' does not exist\n" %(args.input))
		sys.exit()
	else:
		args.input=os.path.realpath(args.input)

	if not args.skipMapping:
		# if fasta provided, the name would be the string before "_R1 or _R2"
		files = glob(args.input+'/*.f*q.gz')

		dnaids=[]
		dnaidsR1=defaultdict(int)
		dnaidsR2=defaultdict(int)

		dnaids = set(sorted([re.split("_R2|_R1", os.path.basename(i))[0] for i in files]))
		for w in [re.split("_R2|_R1", os.path.basename(i))[0] for i in files if '_R1' in i]: dnaidsR1[w] +=1
		for w in [re.split("_R2|_R1", os.path.basename(i))[0] for i in files if '_R2' in i]: dnaidsR2[w] +=1

		if len(dnaids)==0:
			sys.stderr.write("ERROR: Fastq files not found in '%s'\n" %(args.input))
			sys.exit()

		for dnaid in dnaids:
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
		# if bam files provided, the name would be the string before "_"
		dnaids = [(os.path.basename(i)).replace(".bam", "")  for i in glob(args.input+'/*.bam')]
		bai = [(os.path.basename(i)).replace(".bai", "")  for i in glob(args.input+'/*.bai')]
		if len(dnaids) != len(set(dnaids)):
			sys.stderr.write("ERROR: more than one bam file per ID in folder '%s'\n" %(args.input))
			sys.exit()
		elif len(dnaids) != len(bai) or len(dnaids)==0:
			sys.stderr.write("ERROR: check input alignment directory '%s'. Different number of bam and bai files or not existing\n" %(args.input))
			sys.exit()


	# Summary of selected samples

	if args.samples != "all":
		if len(list(set(file_samples).intersection(dnaids))) != len(file_samples):
			sys.stdout.write("\nWARNING - This samples could not be found inside the input directory: \n%s\n\n" %(",".join(set(file_samples) - set(dnaids))))
		if len(set(dnaids))==0:
			sys.stderr.write("ERROR: no fastq.gz/bam files to analyse. \n")
			sys.exit()
		else:
			dnaids = list(set(file_samples).intersection(dnaids))



	sys.stdout.write("\n..................................\n")
	sys.stdout.write("Samples that will be analysed: \n%s" %(",".join(dnaids)))
	sys.stdout.write("\n..................................\n")


	if cat:
		if args.local:
			inputDir = args.output + '/tmp_joinedFastq/'
		else:	
			inputDir = '/scratch/' + os.environ["USER"] + '/' + run + '/'
		if not os.path.exists(inputDir):
			os.makedirs(inputDir)
	else: 
		inputDir = args.input



	#*********** 1. Analysing individual examples if mapping or SNV are specified


	sys.stdout.write("\n......................................................\n")
	sys.stdout.write("  RUNNING MAPPING OR/AND SNV CALLING FOR INDIVIDUAL SAMPLES\n")
	sys.stdout.write("......................................................\n")

	depJobs_list = []

	for sample_name in dnaids: 


			sys.stdout.write("\nAnalysing individual sample '%s' (%s)\n\n" %(sample_name, args.analysis))
			jobid_pr=""

			if args.skipMapping==False:
				
				sampleAnalysis = "mapping"

				#sbatch mapping (n proc)
				myargs_map = [pipelinesPath+"pipeline1_downloadMapping.sh", args.input, args.output, sample_name, str(args.threads), run, "False", str(cat), inputDir,  args.genome, str(args.local), "False", softwarePath]
				jobname_map = "mapping_"+sample_name
				jobid_map=sbatch(jobname_map, args.output, ' '.join(myargs_map), time=args.time, mem=4, threads=args.threads, mail=args.mail, dep='') # mapping alone: just to take advantage of threads
				sys.stdout.write("JOB ID MAPPING: %s\n" %(jobid_map))
				depJobs_list.append(jobid_map)

				# sbatch processing (2 proc - fixed)
				myargs_pr = [tasksPath+"BAMpreprocessing.sh",  str(args.local), run, args.output, sample_name, "False", args.genome]
				#myargs_pr = [pipelinesPath+"pipeline2_SNVprocessingCallingFiltering.sh", args.input, args.output, sample_name, str(args.threads), run, "False", "False", str(cat), inputDir, sampleAnalysis, str(args.cvcf), str(args.skipMapping), args.genome, str(args.local), str(args.pathology), "False", "True", str(removebam), str(args.genefilter), "False", softwarePath]
				print(myargs_pr)
				jobname_pr = "preprocessing_"+sample_name
				jobid_pr=sbatch(jobname_pr, args.output, ' '.join(myargs_pr), time=args.time, mem=args.memory, threads=2, mail=args.mail, dep=jobid_map) # bam processing 
				sys.stdout.write("JOB ID PREPROCESSING: %s\n" %(jobid_pr))
				depJobs_list.append(jobid_pr)



			if args.analysis!="mapping":


			# The haplotypecaller-gvcf-gatk4 workflow runs the HaplotypeCaller tool
			# from GATK4 in GVCF mode on a single sample according to GATK Best Practices.
			# When executed the workflow scatters the HaplotypeCaller tool over a sample
			# using an intervals list file. The output file produced will be a
			# single gvcf file which can be used by the joint-discovery workflow.
			
			# Requirements/expectations :
			# - One analysis-ready BAM file for a single sample (as identified in RG:SM)
			# - Set of variant calling intervals lists for the scatter, provided in a file

				sampleAnalysis = "snv" 
				jobid_hc_list = []
				jobid_annot_list = []

				args.cvcf=True  # haplotype caller in GVCF mode
				args.intervals=True # haplotype caller will be called by intervals
				args.skipMapping == True


				for element in intervals_files:

					print(element)
					args.panel = element

					interval_n = element.split("_temp_")[1].split("_")[0]
					print(interval_n)
					job_name=sampleAnalysis+"-"+interval_n+"_"+sample_name

					output_interval = args.output+"/interval_"+interval_n

					if not os.path.exists(output_interval):
						os.makedirs(output_interval)

					myargs = [args.input, output_interval, sample_name, str(args.threads), run, args.panel, "False", str(cat), inputDir, sampleAnalysis, str(args.cvcf), str(args.skipMapping), args.genome, str(args.local), str(args.pathology), str(args.intervals), "True", str(removebam), str(args.genefilter), "0", "False", str(args.memory*4), softwarePath, interval_n]

					myargs_pipe = [pipelinesPath+"pipeline2_SNVprocessingCallingFiltering.sh"]+myargs # Just calls to SNV processing and calling

					jobid_hc=sbatch(job_name, args.output, ' '.join(myargs_pipe), time=args.time, mem=args.memory, threads=4, mail=args.mail, dep=jobid_pr) 
					jobid_hc_list.append(jobid_hc)
					sys.stdout.write("JOB ID haplotypecaller for interval %s: %s\n" %(interval_n, jobid_hc))
					depJobs_list.append(jobid_hc)


			
				# Merge by interval GVCFs, genotyping and filtering

				myargs = [args.output, sample_name, args.genome, run, str(args.local), str(args.cvcf), softwarePath]
				myargs_pipe = [pipelinesPath+"pipeline4_WGSgenotypingFiltering.sh"]+myargs 

			 	job_name="genotyping"+"_"+sample_name
				depJobs = ':'.join(jobid_hc_list)  	
			 	jobid_genotyp=sbatch(job_name, args.output, ' '.join(myargs_pipe), time=args.time, mem=args.memory, threads=2, mail=args.mail, dep=depJobs) 
			 	sys.stdout.write("JOB ID GENOTYPING: %s\n" %(jobid_genotyp))
			 	depJobs_list.append(jobid_genotyp)




				# Annotation, LOH and processing of output 

				#for chrome in ["chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX"]:
				for chrome	in ["chrX"]:

					myargs = [args.output, sample_name, str(args.threads), run, str(args.local), softwarePath, chrome]
					myargs_pipe = [pipelinesPath+"pipeline5_WGSAnnotation.sh"]+myargs # Just calls to VEP annotation for a interval

				 	job_name="annotation"+chrome+"_"+sample_name
				 	jobid=sbatch(job_name, args.output, ' '.join(myargs_pipe), time=args.time, mem=args.memory, threads=args.threads, mail=args.mail, dep=jobid_genotyp) 
				 	sys.stdout.write("JOB ID SNV ANNOTATION FOR INTERVAL %s: %s\n" %(chrome, jobid))
				 	depJobs_list.append(jobid)
				 	jobid_annot_list.append(jobid)



				# Merge annotated VCFs generated per-interval for individual sample

				myargs = [args.output, sample_name, str(args.threads), run, str(args.local), softwarePath, str(args.genefilter)]
				myargs_pipe = [pipelinesPath+"pipeline4_WGSmergeAnnotationProcessing.sh"]+myargs # Just calls to VEP annotation for a interval

			 	job_name="annotationMergeProc"+"_"+sample_name
				depJobs = ':'.join(jobid_annot_list)  	
			 	jobid=sbatch(job_name, args.output, ' '.join(myargs_pipe), time=args.time, mem=args.memory, threads=2, mail=args.mail, dep=depJobs) 
			 	sys.stdout.write("JOB ID SNV ANNOTATION PROCESSING: %s\n" %(jobid))
			 	depJobs_list.append(jobid)

				

	

	#*********** 2. Remove empty directories


	myargs_remove = [tasksPath+"removeDirs.sh", args.output, str(removebam), bamF]


	sys.stdout.write("\n.................................\n")
	sys.stdout.write("  REMOVING TMP FOLDERS AND FILES \n")
	sys.stdout.write("...................................\n\n")

	if args.local:	
		subprocess.call(myargs_remove, stdout= stdout_f, stderr = stderr_f)

	else:
		depJobs = ':'.join(depJobs_list)
		job_name = "removeDirs_"+run
		jobidREMOVE=sbatch(job_name, args.output, ' '.join(myargs_remove), time=args.time, mem=1, threads=1, mail=args.mail, dep=depJobs, typedep="any")
		depJobs_list.append(jobidREMOVE)
		sys.stdout.write("JOB ID: %s\n" %(jobidREMOVE))
		sys.stdout.write("DEPENDENT JOBS: %s\n" %(depJobs))







	#*********** 3. Running summary



	if args.local:	


		print("................................")
		elapsed = time.time() - fjd_start_time
		fjdtime = str(datetime.timedelta(seconds=elapsed))
		print("TOTAL RUNNING TIME: %s" % fjdtime)
		print("................................")



	else:

		# create a new job which analysed what happend with all the jobs ran during these task

		sys.stdout.write("\n...............................\n")
		sys.stdout.write("  SLURM JOB STATUS SUMMARY \n")
		sys.stdout.write(".................................\n\n")

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



