####################
### SLURM REPORT ###
####################


jobs=$1
outputfile=$2
errfile=$2
outputdir=$3

for job in jobs; do

	# job status
	echo $job"\t" > $output_file
	jobname=sacct -j $job --format="JobName%50" --noheader
	echo $jobname"\t" > output_file
	status=sacct -j $job --format="State" --noheader
	echo $status"\t" > output_file
	if [ “$status != “COMPLETE” ]
		nComplete+=1
		if [ “$status == "FAILED" ]; then
			nComplete+=1
			grep ERROR $3/$jobname.out > $errfile
		fi
		if [ “$status == "CANCELLED" ]; then
			nCancelled+=1
		fi

done

if 
