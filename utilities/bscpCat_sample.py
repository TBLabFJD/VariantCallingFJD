
import sys
import os
from glob import glob
import subprocess
import argparse
import datetime
import shutil
import platform


def joinFastq(basespaceF, fastqFolder, inputF, dnaid, user):

	
	#sys.stdout.write("\nANALYSING SAMPLE %s\n" %(sample))
	if basespaceF=="True":
		dnaid2 = dnaid
		dnaid=dnaid[0:7]

		sampleFolder = fastqFolder+'/'+dnaid2
		os.mkdir(sampleFolder) # check if existing and delete content.
		sys.stdout.write("Downloading fastq files"+"\n")
		
		if user=="False":
			
			config = "."

		else:

			config = user


		cmd = "/home/proyectos/bioinfo/software/bs-cp --io-timeout=60 -q -s //"+ config +"/Projects/"+inputF+"/Samples/"+dnaid2+" "+sampleFolder

		sys.stderr.write(cmd+"\n")
		proc = subprocess.Popen(cmd, shell=True, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
		stdoutTAG, stderrTAG = proc.communicate()
		stdoutTAG=str(stdoutTAG)
		stderrTAG=str(stderrTAG)
		exitcode = proc.returncode
		
		if exitcode != 0:
			sys.stderr.write("ERROR: %s\n" %(stderrTAG))

			if "More than one matching identifier found for" in stderrTAG:
				id = stderrTAG.split('(')[1].split(')')[0]

				#id = stderrTAG.split("\n")[1].split('(',1)[1].split(')')[0]
				sys.stderr.write("WARNING: Downloading sample '%s'\n\n" %(stderrTAG.split('(')[1].split(')')[0].strip()))
				cmd="/home/proyectos/bioinfo/software/bs-cp --io-timeout=60 -q -s //"+ config +"/Projects/"+inputF+"/Samples/"+id+" "+sampleFolder
				exitcode = subprocess.call(cmd, shell=True)
				if exitcode != 0:
					sys.stderr.write("ERROR: 'Something was wrong when downloading sample %s from basespace'.\n" %(stderrTAG.split('(')[1].split(')')[0].strip()))
					sys.exit(1)
			else:
				sys.exit(1)


		foward = sorted(glob(sampleFolder+ '/*_R1*.f*q.gz'))
		reverse = sorted(glob(sampleFolder+ '/*_R2*.f*q.gz'))

	else:
		foward = sorted(glob(inputF+ '/' + dnaid + '*_R1*.f*q.gz'))
		reverse = sorted(glob(inputF+ '/' + dnaid + '*_R2*.f*q.gz'))
	

	if len(foward)==len(reverse) and len(foward)!= 0:
		
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
		sys.exit()


	if basespaceF=="True":
		shutil.rmtree(sampleFolder)

	# if basespaceF=="True":
	# 	shutil.rmtree(sampleFolder)


if __name__ == "__main__":

	import sys
	joinFastq(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5])

