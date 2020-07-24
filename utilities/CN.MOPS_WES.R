### CNV analysis by CN.MOPS
### Author: Lorena de la Fuente
### Date: 07/07/20


rm(list=ls())

#start clock
ptm <- proc.time()
print(ptm)


#**********************#
# package requirements #
#**********************#

library(cn.mops)
library(optparse)


#**************#
#   Arguments  #
#**************#

option_list=list(
  make_option(c("-d","--dir"),type="character", help="Directory with input bam files to analyse."),
  make_option(c('-o','--outputdir'),type="character", help="Output directory."),
  make_option(c('-n','--name'),type="character", help="Project name."),
  make_option(c('-b','--bed'),type="character",default="genome", help="Bed file with genes found in panel"))

opt_parser=OptionParser(option_list = option_list)
opt=parse_args(opt_parser) #list of the args

bamdir <- opt$dir
output_dir <- paste(opt$outputdir,"/copy_number_variation_data/CNMOPS/",sep="")
bedFile <- opt$bed
projectname <- opt$name

# bamdir <- "/mnt/genetica3/Ionut/MCorton/SureSelect_Glaucoma_03/all_bams_new/Sure_123/some"
# output_dir <- "/mnt/genetica/gonzalo/Codex2/"
# bedFile <- "/mnt/genetica/gonzalo/CoNIFER/muestras_Marta/SureSelect_Glaucoma_2018_real_target_manifest_v2_numbered.bed"
# projectname <- "some"

sink(paste(opt$outputdir,"/software_",opt$name,".txt",sep=""),append=TRUE)
print("R SESSION INFO (CODEX2):")
sessionInfo()
sink()

dir.create(output_dir)

setwd(output_dir)



#### Getting the input data from BAM files (also see Section 3.2 and Section 3).


bamFiles <- list.files(bamdir, pattern = '*.bam$', full.names=T)

mysamplenames <- sapply(basename(bamFiles), function(x) sub("_alignment.bam", "",strsplit(x, split = "_",  fixed = T)[[1]][1]))

segments <- read.table(bedFile, sep="\t",as.is=TRUE)
gr <- GRanges(segments[,1],IRanges(segments[,2],segments[,3]))



X <- getSegmentReadCountsFromBAM(bamFiles,GR=gr, sampleNames=mysamplenames)





#### Running the algorithm (also see Section 4.2).

resCNMOPS <- exomecn.mops(X)
resCNMOPS <- calcIntegerCopyNumbers(resCNMOPS)



save(resCNMOPS, file = paste('CNMOPS',"_results.RData", sep=''))



segm <- as.data.frame(segmentation(resCNMOPS))
CNVs <- as.data.frame(cnvs(resCNMOPS))
CNVRegions <- as.data.frame(cnvr(resCNMOPS))

write.table(segm,file="segmentation.csv",row.names = FALSE, sep="\t", quote=F)
write.table(CNVs,file="CNMOPS_results.txt",row.names = FALSE, sep="\t", quote=F)
write.table(CNVRegions,file="cnvr.csv",row.names = FALSE, sep="\t", quote=F)



