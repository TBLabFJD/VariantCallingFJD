####################
### SLURM REPORT ###
####################


#jobs=$1

jobs="5516866 5516867 5516868 5516869 5516870 5516871"
#outputfile=$2
#errfile=$2
#outputdir=$3


for job in jobs; do

	# job status
	echo $job"\t" 
	jobname=sacct -j $job --format="JobName%50" --noheader
	echo $jobname"\t" 
	status=sacct -j $job --format="State" --noheader
	echo $status"\t" 
	# if [ “$status != “COMPLETE” ]
	# 	nComplete+=1
	# 	if [ “$status == "FAILED" ]; then
	# 		nComplete+=1
	# 		grep ERROR $3/$jobname.out 
	# 	fi
	# 	if [ “$status == "CANCELLED" ]; then
	# 		nCancelled+=1
	# 	fi
	# fi
done
 
