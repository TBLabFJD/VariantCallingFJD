#!/bin/sh

##################################################
##  FJD pipeline - Copy Number Variants Calling ##
##################################################


# Arguments
MDAP=$1
HCGVCFD=$2
CNV="${MDAP}/copy_number_variation_data"
SAMPLEFILE=$3
run=$4
threads=$5
panel=$6
window=$7
utilitiesPath=$8
local=$9
methods=${10}
depth=${11}
genefilter=${12}
genome=${13}
# analysis=${14}


# REQUIRED MODULES:


if [ "$local" != "True" ]; then

	module load bedtools/2.27.0
	module load R/R
	module load samtools
	module load vep/release98
	module load annotsv/2.2


	fai="/home/proyectos/bioinfo/references/hg19/ucsc.hg19_convading.fasta.fai"
	VEP="/usr/local/bio/vep/vep"
 	FILTER_VEP='/usr/local/bio/vep/filter_vep'
	VEP_CACHE='/usr/local/bio/vep/t/testdata/cache/homo_sapiens'
	VEP_FASTA="/home/proyectos/bioinfo/references/VEPfasta/Homo_sapiens.GRCh37.dna.primary_assembly.fa"
	PLUGIN_DIR=/usr/local/bio/vep/plugins-95/VEP_plugins-release-95
	PLUGIN_DBS="/home/proyectos/bioinfo/references/VEPdbs"
	dbNSFP_DB="${PLUGIN_DBS}/dbNSFP3.5a_hg19.gz"

	softwareFile="${MDAP}/software_${run}.txt"
	echo "COPY NUMBER VARIANT CALLING:" >> ${softwareFile}
	module list 2>> ${softwareFile}

	export ANNOTSV="/usr/local/bioinfo/annotsv/2.2"
	alias annotsv=$ANNOTSV/bin/AnnotSV/AnnotSV.tcl

else

	VEP_CACHE='/mnt/genetica3/marius/pipeline_practicas_marius/software/variant_effect_predictor/.vep'
	VEP_FASTA="${VEP_CACHE}/homo_sapiens/93_GRCh37/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa"
	VEP="/mnt/genetica3/marius/pipeline_practicas_marius/software/variant_effect_predictor/ensembl-vep/vep"
	fai="/mnt/genetica3/marius/pipeline_practicas_marius/hg19bundle/ucsc.hg19_convading.fasta.fai"

	softwareFile="${MDAP}/software_${run}.txt"
	echo "COPY NUMBER VARIANT CALLING:" >> ${softwareFile}
	printf "\nVEP VERSION\n" >> ${softwareFile}
	$VEP --help | grep "Versions:" -A 5 | tail -n4 >> ${softwareFile}

fi



echo "mkdir ${MDAP}/copy_number_variation_data"
mkdir $CNV


#####################
## QUALITY CONTROL ##
#####################

# sort, remove small regions and create new window size

if echo "$methods" | grep -q "QC"; then


	if [ "$panel" != "genome" ]; then

		printf "\n\nSorting and removing small targets in bed file\n"
		panel_10bp=${MDAP}/$(basename "$panel" .bed)_${run}_10bp.bed
		echo -e "\n"

		echo `wc -l $panel`
		awk '{if(($3-$2)>10 && $1!="chrX" && $1!="chrY" && $1!="Y" && $1!="X"){print $0}}' $panel | sort -V -k1,1 -k2,2 > $panel_10bp
		#awk '{if(($3-$2)>10){print $0}}' $panel | sort -V -k1,1 -k2,2 > $panel_10bp

		echo `wc -l $panel`
		echo `wc -l $panel_10bp`
		printf "\nDONE"
		
		if [ "$window" != "False" ]; then
			printf "\n\nGenerating sliding windows across the BED file"
			panel_out=${MDAP}/$(basename "$panel_10bp" .bed)_125window.bed
			python $utilitiesPath/CNV_windowSize.py $panel_10bp 125 ${panel_out}_unsorted 
			sort -V -k1,1 -k2,2 ${panel_out}_unsorted | uniq > $panel_out
			rm ${panel_out}_unsorted 
		fi
		

		if [ "$depth" != "True" ]; then

			# remove low coverage samples

			bamfiles=$CNV/bamfiles_${run}.txt
			randompanel=${CNV}/$(basename "$panel_10bp" .bed)_covtmp
			samplecov=${CNV}/sample_cov_${run}.txt

			start=`date +%s`
			printf "\n\nComputing genome coverage"
			a=`wc -l $panel_10bp | cut -d " " -f 1`
			N=$((a / 200))
			shuf -n $N $panel_10bp > $randompanel
			ls $HCGVCFD/*bam > $bamfiles
			s1="$?"
			samtools depth -f $bamfiles -b $randompanel > $samplecov
			
			if  [ "$s1" = "0"  ] ; then
				printf '\nEXIT STATUS: 0'
				printf '\nComputing genome coverage for '${run}' run DONE' 
			else
				printf "ERROR: PROBLEMS WITH GENOME COVERAGE ESTIMATION"
				exit 1
			fi

			printf "\n\nFiltering samples according to coverage"

			Rscript $utilitiesPath/regionCoverage.R -f $samplecov -o $CNV -b $bamfiles
			s2="$?"
			
			if  [ "$s2" = "0"  ] ; then
				printf '\nEXIT STATUS: 0'
				printf '\nFiltering samples for '${run}' run DONE' 
				end=`date +%s`
				runtime=$((end-start))
				printf '\nCoverage filtering executing time: '$runtime 
				rm $bamfiles
				rm $randompanel
				rm $samplecov
			else
				printf "ERROR: PROBLEMS WITH SAMPLE FILTERING"
				exit 1
			fi

		fi 
	fi
	

fi





#####################
## COMPUTING CNVs ###
#####################


if [ "$panel" != "genome" ]; then
	if [ "$window" != "False" ]; then
		finalpanel=${MDAP}/$(basename "$panel" .bed)_${run}_10bp_125window.bed
	else
		finalpanel=${MDAP}/$(basename "$panel" .bed)_${run}_10bp.bed
	fi
fi




if echo "$methods" | grep -q "ED"; then
	start=`date +%s`
	printf "\nRunning ExomeDepth"
	#echo Rscript $utilitiesPath/exomeDepth.R -d $HCGVCFD -o $MDAP -b $finalpanel -n $run -s $SAMPLEFILE
	Rscript $utilitiesPath/exomeDepth.R -d $HCGVCFD -o $MDAP -b $finalpanel -n $run -s $SAMPLEFILE

	if [ "$?" = "0" ]; then
		printf '\nEXIT STATUS: 0'
		printf '\nExomeDepth CNV calling for '${run}' DONE'

	else
		printf "\nERROR: PROBLEMS WITH ExomeDepth CNV calling"
		exit 1
	fi

	end=`date +%s`
	runtime=$((end-start))
	echo "ExomeDepth took $(($runtime / 60)) minutes and $(($runtime % 60)) seconds.\n"
fi




if echo "$methods" | grep -q "CN"; then
	start=`date +%s`
	printf "\nRunning CoNVaDING"
	python $utilitiesPath/CoNVading_pipeline.py $HCGVCFD $finalpanel $MDAP $run $utilitiesPath $fai

	if [ "$?" = "0" ]; then
		printf '\nEXIT STATUS: 0'
		printf '\nCoNVading CNV calling for '${run}' DONE'

	else
		printf "\nERROR: PROBLEMS WITH CoNVading CNV calling"
		exit 1
	fi

	end=`date +%s`
	runtime=$((end-start))
	echo "CoNVaDING took $(($runtime / 60)) minutes and $(($runtime % 60)) seconds.\n"
fi




if echo "$methods" | grep -q "C2"; then
	start=`date +%s`
	printf "\nRunning CODEX2"
	Rscript $utilitiesPath/CODEX2.R -d $HCGVCFD -o $MDAP -b $finalpanel -n $run

	if [ "$?" = "0" ]; then
		printf '\nEXIT STATUS: 0'
		printf '\nCODEX2 CNV calling for '${run}' DONE'

	else
		printf "\nERROR: PROBLEMS WITH CODEX2 CNV calling"
		exit 1
	fi

	end=`date +%s`
	runtime=$((end-start))
	echo "CODEX2 took $(($runtime / 60)) minutes and $(($runtime % 60)) seconds.\n"
fi




if echo "$methods" | grep -q "PM"; then
	start=`date +%s`
	printf "\nRunning panelcn.MOPS"
	Rscript $utilitiesPath/panelcnMops.R -d $HCGVCFD -o $MDAP -b $finalpanel -n $run
	end=`date +%s`
	runtime=$((end-start))
	echo "panelcn.MOPS took $(($runtime / 60)) minutes and $(($runtime % 60)) seconds.\n"
fi


###########################################
## DOWNSTREAM ANALYSIS OF DETECTED CNVs ###
###########################################


if echo "$methods" | grep -q "MA"; then

	start=`date +%s`

	printf "\nRunning CNV_result_mixer and gene filtering (optional)"
	Rscript $utilitiesPath/CNV_result_mixer.R -o $MDAP -n $run -b $finalpanel -s $SAMPLEFILE -d $HCGVCFD -f $genefilter
	

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
	

	start=`date +%s`

	printf "\nTabulated CNV results to VCF format"
	Rscript $utilitiesPath/CNV_tsv2vcf.sh $MDAP $run $genome

	end=`date +%s`
	runtime=$((end-start))
	printf '\nExecuting time: '$runtime 



	### IMPORTANT
	### BACK BAM FILES TO ORIGINAL POSITION
	



	#####################
	## CNVs ANNOTATION ##
	#####################

	printf '\n\nCNV annotation... '

	annotsv -SVinputFile  $CNV/${run}"_combined.vcf" -outputfile $CNV/${run}"_combinedAnnotated.tsv"



	if [ "$?" = "0" ]; then
		printf '\nEXIT STATUS: 0'
		printf "\nSV ANNOTATION FOR FILE $CNV/${run}"_combined.txt" DONE"

	else
		printf "\nERROR: PROBLEMS WITH ANNOTSV ANNOTATION"
		exit 1
	fi





	# ##############################
	# ## FINAL CNV ANNOTATED FILE ##
	# ##############################


	echo -e "\nMerging exomeDepth and VEP output files..."


	Rscript $utilitiesPath/CNVfinalMerge.R -d $MDAP -r $run


	if [ "$?" = "0" ]; then
		printf '\nEXIT STATUS: 0'
		printf "\nSV ANNOTATION FOR FILE $CNV/${run}"_combinedAnnotated.tsv" DONE"

	else
		printf "\nERROR: PROBLEMS WITH MERGING MIXER AND ANNOTATED FILE"
		exit 1
	fi



fi


