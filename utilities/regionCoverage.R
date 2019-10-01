### Filtering low coverage samples for downstream CNV analysis
### Author: Lorena de la Fuente
### Date: 14/01/19


library(optparse)


# arguments

option_list=list(
  make_option(c("-f","--file"),type="character", help="Coverage file."),
  make_option(c("-b","--bamfile"),type="character", help="Names bam files."),
  make_option(c('-o','--outputdir'),type="character", help="Output directory."))

opt_parser=OptionParser(option_list = option_list)
opt=parse_args(opt_parser) #list of the args

resultsFile = paste(opt$outputdir, "/minCovFilterResults.txt", sep = "")


# coverage
df = NULL

bams = read.delim(opt$bamfile, header=F, stringsAsFactors = F)[,1]
coverage = read.delim(opt$file, header = F, stringsAsFactors = F)
coverage = coverage[,3:ncol(coverage)]
colnames(coverage) = bams

bamdir = strsplit(opt$file, split = "/",  fixed = T)[[1]]
bamdir = paste(bamdir[1:length(bamdir)-1] , collapse="/")

for (n in 1:ncol(coverage)){
  #cmd = paste("bedtools coverage -b",file,"-a",opt$panel, "-d", sep=" ")
  #cat('\nBedtools coverage for sample', file, sep = " ")
  #cov = system(cmd, intern = TRUE)
  file = colnames(coverage)[n]
  filewopath = tail(strsplit(file, split = "/",  fixed = T)[[1]], n=1)
  sample =  sub(".bam", "",strsplit(filewopath, split = "_",  fixed = T)[[1]][1])
  cov = coverage[,n]
  s = summary (as.numeric (cov))
  m = median (as.numeric (cov))
  if(m>20){Filter="PASSED"}else{Filter="FAILED"}
  df = rbind(df,c(sample, as.character(s), as.character(Filter), file))
  
}

colnames(df) = c("Sample", names(s), "Min.Cov.Test", "File")

write.table(df, file=resultsFile, append = T, row.names = F, col.names = T, quote = F, sep="\t")


if (any(data.frame(df)$Min.Cov.Test=="FAILED")){
  setwd(bamdir)
  dir.create("lowCov_bams")
  filtfiles = df[data.frame(df)$Min.Cov.Test=="FAILED","File"]
  for (f in filtfiles){
    st=paste(sub('\\.bam$', '', f),"*", sep="")
    system(paste("mv",st,"lowCov_bams/", sep=" "))
    }
  }


cat('\nSee', resultsFile, 'to check filtered low-coverage samples', '\n', sep = " ")

