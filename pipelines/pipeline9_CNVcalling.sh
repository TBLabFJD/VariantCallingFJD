#!/bin/sh

##################################################
##  Copy Number Variants Calling and Annotaton ##
##################################################


# Arguments
INPUT=$1
MDAP=$2
skipmapping=$3
CNV="${MDAP}/cnvs"
mkdir $CNV
samplenames=$4
run=$5
threads=$6
panel=$7
window=$8
softwarePath=$9
local=${10}
methods=${11}
genefilter=${12}
genome=${13}
mycountthreshold=${14}
sexchr=${15}
tasksPath=${softwarePath}/tasks




# REQUIRED MODULES:


#if [ "$local" != "True" ]; then

module load bedtools/2.27.0
module load R/R
module load samtools/1.9
module load annotsv/2.2
module load perl/5.28.0
source ${softwarePath}/pipeline.config

fai=${fai_path}
#fai="/home/proyectos/bioinfo/references/hg19/ucsc.hg19_convading.fasta.fai"

softwareFile="${MDAP}/software_${run}.txt"
echo "COPY NUMBER VARIANT CALLING:" >> ${softwareFile}
module list 2>> ${softwareFile}
export ANNOTSV=${annotsv_path}
# export ANNOTSV="/usr/local/bioinfo/annotsv/2.2"

# else

# 	fai="/mnt/genetica3/marius/pipeline_practicas_marius/hg19bundle/ucsc.hg19_convading.fasta.fai"

# 	softwareFile="${MDAP}/software_${run}.txt"
# 	echo "COPY NUMBER VARIANT CALLING:" >> ${softwareFile}
# 	printf "\nVEP VERSION\n" >> ${softwareFile}
# 	$VEP --help | grep "Versions:" -A 5 | tail -n4 >> ${softwareFile}

# fi






#####################
## SET BAM FOLDER  ##
#####################


# Define bam folder

if [ "$skipmapping" != "True" ]; then bamF=$MDAP/bams; else bamF=${INPUT}; fi






#####################
## QUALITY CONTROL ##
#####################

# calculate coverage, remove small regions and create new window  size if required.

if echo "$methods" | grep -q "QC"; then


	# Check bed file has four columns
	printf "\n\nChecking bed file has 4 columns..."

	n=$(awk '{print NF}' $panel | sort -nu | wc -l) 
	ncol=$(awk '{print NF}' $panel | sort -nu | tail -n1) 
	if [ "$n" != "1" ] || [ "$ncol" != "4" ]; then echo "ERROR: Target file (panel file) has not 4 columns. Check bed file format."; exit 1; else printf "\nDONE"; fi


	# Sorting and removing small targets - bed file
	printf "\n\nSorting and removing small targets in bed file..."
	panel_10bp=${MDAP}/$(basename "$panel" .bed)_${run}_10bp.bed
	
	if [ "$sexchr" != "False" ]; then
		awk '{if(($3-$2)>10 && ($1=="chrX" || $1=="X") ){print $0}}' $panel | sort -V -k1,1 -k2,2 > $panel_10bp
		#awk '{if(($3-$2)>10 && ($1=="chrX" || $1=="chrY" || $1=="Y" || $1=="X") ){print $0}}' $panel | sort -V -k1,1 -k2,2 > $panel_10bp
	else
		awk '{if(($3-$2)>10 && $1!="chrX" && $1!="chrY" && $1!="Y" && $1!="X"){print $0}}' $panel | sort -V -k1,1 -k2,2 > $panel_10bp

	fi

	if [ "$?" = "1" ]  ; then echo "ERROR: Check targets in panel file."; exit 1;  else printf "\nDONE"; fi
	

	# Divide targets in windowns
	if [ "$window" != "False" ]; then
		printf "\n\nGenerating sliding windows across the BED file..."
		panel_out=${MDAP}/$(basename "$panel_10bp" .bed)_125window.bed
		python $tasksPath/CNV_windowSize.py $panel_10bp 125 ${panel_out}_unsorted 
		if [ "$?" = "1" ]  ; then echo "ERROR: Problem when generating sliding window across targets."; exit 1;  else printf "\nDONE"; fi
	
		sort -V -k1,1 -k2,2 ${panel_out}_unsorted | uniq > $panel_out
		rm ${panel_out}_unsorted 
	fi
	

	#Compute coverage level by sample 

	printf "\n\nChecking if test samples passed the coverage QC:\n"
	for bam in $bamF/*bam; do

		sample="$(basename $bam .bam)"
		$tasksPath/mosdepth.sh $local $run $MDAP $sample $bam $panel 0 cnv $mycountthreshold

	done

	printf "\n\nGenerating coverage plot..\n"
	python $tasksPath/plot-dist.py ${MDAP}/qc/*_qc_cnv.mosdepth.region.dist.txt -o ${MDAP}/qc/mosdepth.region.dist_cnv_${run}.html
	if [ "$?" = "1" ]  ; then echo "ERROR: Problems generating Plot Distribution."; exit 1;  else printf "\nDONE"; fi
	rm ${MDAP}/qc/*qc_cnv*


	#Create folder of samples passing QC parameters if any was FAILED
	
	n=$(cut -f4 $MDAP/qc/minCovFilterResults_cnv_${run}.txt | uniq -c | wc -l)
	status=$(cut -f4 $MDAP/qc/minCovFilterResults_cnv_${run}.txt | uniq -c | awk '{print $2}')
	if [ "$n" != "1" ] || [ "$status" != "PASSED" ]; then 


		echo "WARNING: Bam files will be analyzed from $CNV/passed_bams folder"
		printf "\n\nSelecting bam files passing coverage QC...\n"
		bamF_pass=$CNV/passed_bams
		mkdir $bamF_pass
		qc_file=$MDAP/qc/minCovFilterResults_cnv_${run}.txt
		samples=$(awk '{if($4=="PASSED"){print $6}}' $qc_file | sed 's/.bam$//g')
		s1="$?"
		for file in $samples; do cp ${file}*bam ${file}*bai $bamF_pass/.; done
		s2="$?"
		if [ "$s1" = "1" ] || [ "$s2" = "1" ] ; then exit 1;  else printf "\nDONE"; fi

	
	else

		printf "\nAll bam files passed the QC."

	fi

	exit 0

fi








#####################
## COMPUTING CNVs ###
#####################




# Defining the panel (generated during QC)

if [ "$panel" != "genome" ]; then
	if [ "$window" != "False" ]; then
		finalpanel=${MDAP}/$(basename "$panel" .bed)_${run}_10bp_125window.bed
	else
		finalpanel=${MDAP}/$(basename "$panel" .bed)_${run}_10bp.bed
	fi
fi


# Defining the bam folder (created during QC)
n=$(cut -f4 $MDAP/qc/minCovFilterResults_cnv_${run}.txt | uniq -c | wc -l)
status=$(cut -f4 $MDAP/qc/minCovFilterResults_cnv_${run}.txt | uniq -c | awk '{print $2}')
if [ "$n" != "1" ] || [ "$status" != "PASSED" ]; then echo "WARNING: Bam files will be analyzed from $CNV/passed_bams folder"; bamF_pass=$CNV/passed_bams; else bamF_pass=$bamF; fi



if echo "$methods" | grep -q "ED"; then
	start=`date +%s`
	printf "\nRunning ExomeDepth"
	Rscript $tasksPath/exomeDepth.R -d $bamF_pass -o $MDAP -b $finalpanel -n $run

	if [ "$?" = "0" ]; then
		printf '\nEXIT STATUS: 0'
		printf '\nExomeDepth CNV calling for '${run}' DONE'

	else
		printf "\nERROR: PROBLEMS WITH ExomeDepth CNV calling"
		exit 1
	fi

	end=`date +%s`
	runtime=$((end-start))
	echo "\nExomeDepth took $(($runtime / 60)) minutes and $(($runtime % 60)) seconds.\n"
fi




if echo "$methods" | grep -q "CN"; then
	start=`date +%s`
	printf "\nRunning CoNVaDING"
	python $tasksPath/CoNVading_pipeline.py $bamF_pass $finalpanel $MDAP $run $tasksPath $fai

	if [ "$?" = "0" ]; then
		printf '\nEXIT STATUS: 0'
		printf '\nCoNVading CNV calling for '${run}' DONE'

	else
		printf "\nERROR: PROBLEMS WITH CoNVading CNV calling"
		exit 1
	fi

	end=`date +%s`
	runtime=$((end-start))
	echo "\nCoNVaDING took $(($runtime / 60)) minutes and $(($runtime % 60)) seconds.\n"
fi




if echo "$methods" | grep -q "C2"; then
	start=`date +%s`
	printf "\nRunning CODEX2"
	Rscript $tasksPath/CODEX2.R -d $bamF_pass -o $MDAP -b $finalpanel -n $run

	if [ "$?" = "0" ]; then
		printf '\nEXIT STATUS: 0'
		printf '\nCODEX2 CNV calling for '${run}' DONE'

	else
		printf "\nERROR: PROBLEMS WITH CODEX2 CNV calling"
		exit 1
	fi

	end=`date +%s`
	runtime=$((end-start))
	echo "\nCODEX2 took $(($runtime / 60)) minutes and $(($runtime % 60)) seconds.\n"
fi




if echo "$methods" | grep -q "PM"; then
	start=`date +%s`
	printf "\nRunning panelcn.MOPS"
	Rscript $tasksPath/panelcnMops.R -d $bamF_pass -o $MDAP -b $finalpanel -n $run
	end=`date +%s`
	runtime=$((end-start))
	echo "\npanelcn.MOPS took $(($runtime / 60)) minutes and $(($runtime % 60)) seconds.\n"
fi








if echo "$methods" | grep -q "MA"; then

	start=`date +%s`



	###########################################################
	## MERGING CNV CALLING RESULTS FROM ALTERNATIVE METHODS ###
	###########################################################


	printf "\nMerging CNV calling results from alternative methods\n"
	Rscript $tasksPath/CNV_result_mixer.R -o $MDAP -n $run -b $finalpanel -s $samplenames -d $bamF 
	

	if [ "$?" = "0" ]; then
		printf '\nEXIT STATUS: 0'
		printf '\nCNV mix for '${run}' DONE'

	else
		printf "\nERROR: PROBLEMS WITH CNV mixer"
		exit 1
	fi

	end=`date +%s`
	runtime=$((end-start))
	printf '\nExecuting time: '$runtime 
	


	################################
	## CNV RESULTS TO VCF FORMAT ###
	################################


	start=`date +%s`

	printf "\n\nTabulated CNV results to VCF format"
	$tasksPath/CNV_tsv2vcf.sh $MDAP $run $genome
	if [ "$?" = "0" ]; then
		printf '\nEXIT STATUS: 0'
		printf '\nCNV FORMAT CHANGE for '${run}' DONE'

	else
		printf "\nERROR: PROBLEMS WITH CNV RESULTS TO VCF"
		exit 1
	fi

	end=`date +%s`
	runtime=$((end-start))
	printf '\nExecuting time: '$runtime 
	


	#####################
	## CNVs ANNOTATION ##
	#####################

	printf '\n\nCNV annotation...\n'

	${annotsv_path}/bin/AnnotSV/AnnotSV.tcl -SVinputFile  $CNV/${run}.combined.vcf -outputfile $CNV/${run}.combinedAnnotated.tsv

	if [ "$?" = "0" ]; then
		printf '\nEXIT STATUS: 0'
		printf "\nSV ANNOTATION FOR FILE "${run}".combined.txt DONE"

	else
		printf "\nERROR: PROBLEMS WITH ANNOTSV ANNOTATION"
		exit 1
	fi





	# ##############################
	# ## FINAL CNV ANNOTATED FILE ##
	# ##############################


	echo -e "\nMerging exomeDepth and VEP output files..."


	Rscript $tasksPath/CNVfinalMerge.R -d $MDAP -r $run


	if [ "$?" = "0" ]; then
		printf '\nEXIT STATUS: 0'
		printf "\nFINAL MERGE FOR FILE "${run}".combinedAnnotated.tsv DONE"

	else
		printf "\nERROR: PROBLEMS WITH MERGING MIXER AND ANNOTATED FILE"
		exit 1
	fi



	####################
	## GENE FILTERING ##
	####################


	echo -e "\nFiltering results by genefilter..."


	if [ "$genefilter" != "False" ]; then

		head -n1 $CNV/${run}.final.txt > $CNV/${run}.final.genelist.txt
		awk -F "\t" 'FNR==NR { a[$1]; next } $3 in a || $18 in a'  $genefilter $CNV/${run}.final.txt >> $CNV/${run}.final.genelist.txt

	fi

	if [ "$?" = "0" ]; then
		printf '\nEXIT STATUS: 0'
		printf "\nGENE FILTERING BY INPUT USER LIST DONE"

	else
		printf "\nERROR: PROBLEMS WITH GENE FILTERING"
		exit 1
	fi



	# Removing temporal folder with passed files

	rm -r $bamF_pass



fi


