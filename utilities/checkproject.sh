#!/bin/bash

####################
### SLURM REPORT ###
####################
# jobs=$1 # all jobs in a variable
# outputfile=$2
# errfile=$2
# outputdir=$3



#jobs=""
#jobid=`sbatch --account=bioinfo_serv --partition=bioinfo sleep.sh 20`
#jobid=($jobid)
#id=${jobid[3]}
#jobs=$jobs" "$id

#jobid=`sbatch --account=bioinfo_serv --partition=bioinfo sleep.sh 60`
#jobid=($jobid)
#id=${jobid[3]}
#jobs=$jobs" "$id

jobs="5350274 5350275"
echo $jobs


# output_file="/home/proyectos/bioinfo/fjd/ionut_lorenaRun/test_check_project_status"
# jobs="5304555 5304556 5304557 5304558 5304559 5304560 5304561 5304562 5304563 5304564 5304565 5304566 5304567 5304568 5304569 5304570 5304571 5304572 5304573 5304574 5304575 5304576 5304577 5304578 5304579 5304580 5304581 5304582 5304583 5304584 5304585 5304586 5304587 5304588 5304589 5304590 5304591 5304592 5304593 5304594 5304595 5304596 5304597 5304598 5304599 5304600 5304601 5304602 5304603"

# jobs in a dictionary with key as type of job.

nFailed=0
nCancelled=0
nComplete=0

for job in $jobs; do

	# job status
	jobname=`sacct -j $job --format="JobName%50" --noheader`
	status=`sacct -j $job --format="State" --noheader`
	echo "hola"
	echo -e $job"\t"${jobname}"\t"${status}
	echo "adios"

	if [ "$status" = "COMPLETE" ]; then
		nComplete=$((nComplete+1))
	else
		echo "hola"
		if [ "$status" = "FAILED" ]; then 
			nFailed=$((nFailed+1))
			grep ERROR $3/${jobname}.out > $output_file
		fi
		if [ "$status" = "CANCELLED" ]; then
			nCancelled=$((nCancelled+1))
		fi
	fi

done


echo $nCancelled
echo $nComplete
echo $nFailed



# echo 'ndatasets npvm nerror nmemory' > summary.txt

# for folder in `ls`; do

#         echo $folder
#         ndatasets=`ls -l $folder/mapping*out | wc -l`;
#         pvm=`grep "POST-VEP MODIFICATIONS" $folder/mapping*out | wc -l`;
#         nerror=`tail -n2 $folder/*err | grep "ERROR" | wc -l`;
#         nmemory=`tail -n2 $folder/*err | grep "CANCELLED" | wc -l`;
#         echo $ndatasets $pvm $nerror $nmemory;


# done >> summary.txt
