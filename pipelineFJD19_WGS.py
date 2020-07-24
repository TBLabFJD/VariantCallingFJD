#!/usr/bin/env python

#######################################################################
# Script: Run FJD analysis for samples stored in basespace or locally
# Author: Lorena de la Fuente 
# Date: 24-01-2019
#######################################################################

import sys
import os
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

sys.path.insert(0, utilitiesPath)


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
	parser = argparse.ArgumentParser(description="FJD genotyping analysis for samples from Basespace")
	parser.add_argument('-i', '--input', help='\t\t Local directory with fastq.gz/bam input files or project name within Basespace', required=True)
	parser.add_argument('-o', '--output', dest="output",help='\t\tOutput directory name.', required=True)
	parser.add_argument('-s', '--samples' , help='\t\tCSV file listing sample names within second column', required=False, default="all")
	parser.add_argument('-n', '--name', help='\t\tName for job. RECOMMENDED', required=False)
	parser.add_argument('-a', '--analysis', help='\t\tType of analysis to run', required=False, choices={"mapping", "snp", "cnv", "all"}, default="snp")
	parser.add_argument('-c', '--cvcf', help='\t\tCombined genotyping. Number of samples must be higher than 2', required=False, action='store_true')
	parser.add_argument('-p', '--panel', help='\t\tBed with panel regions. Mandatory if CNV analysis', required=False)
	parser.add_argument('-l', '--local', help='\t\tRun in local', required=False, action='store_true')
	parser.add_argument('-t', '--threads', help='\t\tNumber of threads', type=int, required=False, default=1)
	parser.add_argument('-M', '--memory', help='\t\tNumber of GBs for VEP running', type=int, required=False, default=5)
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
	parser.add_argument('-D', '--depth', help="Skip filtering of samples during CNV calling by read coverage.", required=False, action='store_true')
	parser.add_argument('-w', '--window', help="125 window size for CNV intervals", required=False, action='store_true')
	parser.add_argument('-B', '--remove_bam', help="Remove bam files", required=False, action='store_true')
	parser.add_argument('-d', '--duplicates', help="Data is not dedupped. By default, duplicates are just marked", required=False, action='store_true')
	parser.add_argument('-u', '--basemountuser', help="Specify alternative basemount user", required=False, default=False)
	parser.add_argument('-b', '--basespace', help='\t\tTake samples from Basespace', required=False, action='store_true')


	fjd_start_time = time.time()
	args = parser.parse_args()


	# label for current lab

	if args.name!=None:
		run = args.name+"_"+now
	else:
		run = now


	# panel 

	args.panel = "genome"


	## defining stdout and stderr file names 

	stdoutAll="FJD_"+run+".out"
	sys.stdout = open(stdoutAll, 'w', 0)
	stderrAll="FJD_"+run+".err"
	sys.stderr = open(stderrAll, 'w', 0)



	## pipeline starts

	sys.stdout.write("\nFJD ANALYSIS FOR WHOLE GENOME SEQUENCING\n")
	sys.stdout.write("\nChecking arguments...\n")



	# checking output dir

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
		
		if len(dnaids) > 1:
			sys.stderr.write("ERROR: more than one sample in folder '%s'. Just one sample must be provided for WGS analysis\n" %(args.input))
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
		dnaids = [(os.path.basename(i)).split("_")[0].replace(".bam", "")  for i in glob(args.input+'/*.bam')]
		bai = [(os.path.basename(i)).split("_")[0].replace(".bai", "")  for i in glob(args.input+'/*.bai')]
		if len(dnaids) != len(set(dnaids)):
			sys.stderr.write("ERROR: more than one bam file per ID in folder '%s'\n" %(args.input))
			sys.exit()
		elif len(dnaids) != len(bai) or len(dnaids)==0:
			sys.stderr.write("ERROR: check input alignment directory '%s'. Different number of bam and bai files or not existing\n" %(args.input))
			sys.exit()
		if len(dnaids) > 1:
			sys.stderr.write("ERROR: more than one alignment file in folder '%s'. Just one sample must be provided for WGS analysis\n" %(args.input))
			sys.exit()


	# Sample summary
	sample_name = dnaids.pop()	#tested samples




	sys.stdout.write("\n..................................\n")
	sys.stdout.write("Sample that will be analysed: \n%s" %(sample_name))
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



	# Analysing individual examples if mapping or SNP are specified


	jobid_list=[]
	jobid_list_snp=[]
	depJobs=''
	jobid = ''


	sys.stdout.write("\n......................................................\n")
	sys.stdout.write("  RUNNING MAPPING OR/AND SNV CALLING FOR INDIVIDUAL SAMPLES\n")
	sys.stdout.write("......................................................\n")

	sys.stdout.write("\nAnalysing individual sample '%s' (%s)\n\n" %(sample_name, args.analysis))


	if args.skipMapping==False:
		
		sampleAnalysis = "mapping"

		myargs_desc = ["INPUT FOLDER", "OUTPUT DIR", "SAMPLE LABEL", "N THREADS", "RUN LABEL", "PANEL BED FILE", "BASESPACE DOWNLOAD", "CONCATENATION", "INPUT DIRECTORY", "ANALYSIS TYPE", "COMBINED VCF", "SKIPPING MAPPING STEP", "GENOME BUNDLE", "LOCAL", "PATHOLOGY", "INTERVALS", "REMOVE DUPLICATES", "REMOVE BAM FILES", "GENE LIST", "BASESPACE USER", "SOFTPATH"]
		myargs = [args.input, args.output, sample_name, str(args.threads), run, args.panel, str(args.basespace), str(cat), inputDir, sampleAnalysis, str(args.cvcf), str(args.skipMapping), args.genome, str(args.local), str(args.pathology), str(args.intervals), str(args.duplicates), str(removebam), str(args.genefilter), str(args.basemountuser) ,softwarePath]
		[sys.stdout.write("%s: %s\n" %(myargs_desc[i], myargs[i])) for i in range(1,len(myargs)-1)]
		
		myargs_pipe1 = [pipelinesPath+"pipeline1_downloadMapping.sh"]+myargs
		jobname_pipe1 =  "mapping"+"_"+sample_name

		myargs_pipe2 = [pipelinesPath+"pipeline2_SNVprocessingCallingFiltering.sh"]+myargs
		jobname_pipe2 =  "processing_"+sampleAnalysis+"_"+sample_name

		# sbatch mapping (n proc)
		tmpjobid=sbatch(jobname_pipe1, args.output, ' '.join(myargs_pipe1), time=args.time, mem=4, threads=args.threads, mail=args.mail, dep='') # mapping alone: just to take advantage of threads
		sys.stdout.write("JOB ID MAPPING: %s\n" %(tmpjobid))

		# sbatch processing (2 proc - fixed)
		jobid=sbatch(jobname_pipe2, args.output, ' '.join(myargs_pipe2), time=args.time, mem=args.memory, threads=2, mail=args.mail, dep=tmpjobid) # bam processing 
		jobid_list.append(jobid)
		sys.stdout.write("JOB ID PROCESSING: %s\n" %(jobid))



	if args.analysis!="mapping":

	## The haplotypecaller-gvcf-gatk4 workflow runs the HaplotypeCaller tool
	## from GATK4 in GVCF mode on a single sample according to GATK Best Practices.
	## When executed the workflow scatters the HaplotypeCaller tool over a sample
	## using an intervals list file. The output file produced will be a
	## single gvcf file which can be used by the joint-discovery workflow.
	##
	## Requirements/expectations :
	## - One analysis-ready BAM file for a single sample (as identified in RG:SM)
	## - Set of variant calling intervals lists for the scatter, provided in a file

		sampleAnalysis = "snp" 

		args.cvcf=True  # haplotype caller in GVCF mode
		args.intervals=True
		args.skipMapping == True

		# for element in intervals_files:

		# 	print(element)
		# 	args.panel = element

		# 	interval_n = element.split("_temp_")[1].split("_")[0]
		# 	print(interval_n)
		# 	job_name=sampleAnalysis+"_"+interval_n

		# 	output_interval = args.output+"/interval_"+interval_n
		# 	os.makedirs(output_interval)

		# 	myargs = [args.input, output_interval, sample_name, str(args.threads), run, args.panel, str(args.basespace), str(cat), inputDir, sampleAnalysis, str(args.cvcf), str(args.skipMapping), args.genome, str(args.local), str(args.pathology), str(args.intervals), str(args.duplicates), str(removebam), str(args.genefilter), str(args.basemountuser), softwarePath, interval_n]

		# 	myargs_pipe = [pipelinesPath+"pipeline2_SNVprocessingCallingFiltering.sh"]+myargs

		# 	jobid=sbatch(job_name, args.output, ' '.join(myargs_pipe), time=args.time, mem=args.memory, threads=2, mail=args.mail, dep=jobid) 
		# 	jobid_list_snp.append(jobid)
		# 	sys.stdout.write("JOB ID HAPLOTYPE CALLER CALLING FOR INTERVAL %s: %s\n" %(interval_n, jobid))


		
		# Merge GVCFs generated per-interval for the same sample

		# filenames = glob(args.output+"/interval_*/genotyping/*.list")


		# gvcffile = args.output+"/"+os.path.basename(filenames[0])

		# with open(gvcffile, 'w') as outfile:
		#     for fname in filenames:
		#         with open(fname) as infile:
		#             for line in infile:
		#                 outfile.write(line)


		# myargs = [args.output, gvcffile, sample_name,  args.genome]

		# myargs_pipe = [tasksPath+"mergeGVCF.sh"]+myargs



		# job_name="mergeGVCF_"+sample_name
		# jobid=sbatch(job_name, args.output, ' '.join(myargs_pipe), time=args.time, mem=args.memory, threads=2, mail=args.mail, dep=jobid) 
		# jobid_list_snp.append(jobid)
		# sys.stdout.write("JOB ID mergeGVCF: %s\n" %(jobid))



		# Genotyping 

		myargs = [str(args.local), run, args.output, sample_name, args.genome]
		myargs_pipe = [tasksPath+"genotyping.sh"]+myargs

		job_name="genotypeGVCF_"+sample_name
		jobid=sbatch(job_name, args.output, ' '.join(myargs_pipe), time=args.time, mem=args.memory, threads=2, mail=args.mail) 
		jobid_list_snp.append(jobid)
		sys.stdout.write("JOB ID genotyping: %s\n" %(jobid))



		# Filtering SNVs 

		myargs = [str(args.local), run, args.output, sample_name, args.genome, "CNN", str(args.cvcf)]
		myargs_pipe = [tasksPath+"SNVfiltering.sh"]+myargs

		job_name="filteringGVCF_"+sample_name
		jobid=sbatch(job_name, args.output, ' '.join(myargs_pipe), time=args.time, mem=args.memory, threads=2, mail=args.mail, dep=jobid) 
		jobid_list_snp.append(jobid)
		sys.stdout.write("JOB ID genotyping: %s\n" %(jobid))



	sys.stdout.close()
	sys.stderr.close()


if __name__ == "__main__":
	main()



