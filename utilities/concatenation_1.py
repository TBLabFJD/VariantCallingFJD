#!/usr/bin/env python


import sys
import os
from glob import glob
import subprocess
import argparse
from collections import Counter
import datetime
import shutil



def joinFastq(basespaceF, fastqFolder, inputF, sampleFile, analysis, output):


	# checking and reading sample names in sample file (if provided)

	file_samples = list()

	if sampleFile != "all":
		samplesFile = open(sampleFile, "r")
		[file_samples.append((line.split(",")[1]).strip()) for line in samplesFile]
		file_samples = set(file_samples)


	if basespaceF=="True":

		### Mount basespace with basemount
		sys.stdout.write("\nMounting Basespace...\n")
		basespaceDir=output+"/tmp_basespace"
		print(basespaceDir)
		if os.path.isdir(basespaceDir): 
			sys.stderr.write("ERROR: Temporary basespace folder '%s' already created. Change output path.\n" %(basespaceDir))
		else:
			os.mkdir(basespaceDir)
			subprocess.call(["basemount",basespaceDir])
			sys.stdout.write("\nMounting Basespace COMPLETE...\n")

		if not os.path.isdir(basespaceDir+'/Projects/'+inputF): 
			sys.stderr.write("ERROR: Project '%s' does not exist in basespace\n" %(inputF))
			sys.stdout.write("\nUnmounting Basespace...\n")
			subprocess.call(["basemount","--unmount", basespaceDir])			
			sys.exit()
		else:
			dnaids =  sorted([f for f in os.listdir(basespaceDir+'/Projects/'+inputF+'/Samples') if not f.startswith('.')])
	
	else: 
		inputF = os.path.realpath(inputF)
		dnaids = sorted(set([(os.path.basename(i)).split("_")[0].replace(".fastq.gz", "")  for i in glob(inputF+'/*.fastq.gz')]))


	print(dnaids)

	# sample concatenation 

	sample_names = list()
	for dnaid in dnaids:
		
		#sys.stdout.write("\nANALYSING SAMPLE %s\n" %(sample))
		if basespaceF=="True":
			dnaid2 = dnaid
			dnaid=dnaid[0:7]
			print(dnaid)

		if  analysis in ["cnv","all"] or sampleFile=="all" or dnaid in file_samples:
			
			if basespaceF=="True":
				print("antes rsync")
				print(datetime.datetime.now())
				os.mkdir(fastqFolder+'/'+dnaid2) # check if existing and delete content.
				args = ["rsync", "--progress", "--chmod=777", "-r", basespaceDir+'/Projects/'+inputF+'/Samples/' + dnaid2 + "/Files", fastqFolder+'/'+dnaid2]
				subprocess.call(args)
				print("despues rsync")
				print(datetime.datetime.now())

				foward = sorted(glob(fastqFolder+"/"+dnaid2+ '/Files/*R1*.fastq.gz'))
				reverse = sorted(glob(fastqFolder+"/"+dnaid2+ '/Files/*R2*.fastq.gz'))

			else:
				foward = sorted(glob(inputF+ '/' + dnaid + '*R1*.fastq.gz'))
				reverse = sorted(glob(inputF+ '/' + dnaid + '*R2*.fastq.gz'))
			
			if len(foward)==len(reverse) and len(foward)!= 0:
				sys.stdout.write("\n- SAMPLE: "+ dnaid+"\n")
				
				### FASTQ CONCATENATION	

				foward_file = fastqFolder + dnaid + '_R1.fastq.gz'
				reverse_file = fastqFolder + dnaid + '_R2.fastq.gz'

				sys.stdout.write(foward_file + " written"+"\n")
				sys.stdout.write(reverse_file +" written"+"\n")

				#sys.stdout.write('Joining forward fastq files of sample ' + dnaid +"\n")
				mycmd = 'cat %s > %s' %(' '.join(foward), foward_file)
				subprocess.call(mycmd, shell=True)

				#sys.stdout.write('Joining reverse fastq files of sample ' + dnaid +"\n")
				mycmd2 = 'cat %s > %s' %(' '.join(reverse), reverse_file)
				subprocess.call(mycmd2, shell=True)

			else:
				sys.stderr.write("ERROR: Not fastq files found for sample '%s' or different number of reverse and foward fastq files.\n" %(dnaid))


			#if basespaceF=="True":
			#	shutil.rmtree(fastqFolder+'/'+dnaid2)


	if basespaceF=="True":

		## Mount basespace with basemount
		sys.stdout.write("\nUnmounting Basespace...\n")
		subprocess.call(["basemount","--unmount", basespaceDir])


if __name__ == "__main__":
    
    import sys
    joinFastq(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4],sys.argv[5], sys.argv[6])


