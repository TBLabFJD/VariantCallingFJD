### CODEX COMPARISON
### Author: Lorena de la Fuente


library(CODEX2)

load("/home/proyectos/bioinfo/lodela/SVbenchmark/CODEX2/copynumber2/copy_number_variation_data/CODEX2/CODEX2_counts_image.RData")

output_dir="/home/proyectos/bioinfo/lodela/SVbenchmark/CODEX2/copynumber2/copy_number_variation_data/CODEX2/"
setwd(output_dir)



##### Quality control 

qcObj <- qc(Y, sampname, ref, cov_thresh = c(20, 10000),
            length_thresh = c(20, 500), mapp_thresh = 0.9,
            gc_thresh = c(20, 80))
Y_qc <- qcObj$Y_qc; sampname_qc <- qcObj$sampname_qc
ref_qc <- qcObj$ref_qc; qcmat <- qcObj$qcmat; gc_qc <- ref_qc$gc
write.table(qcmat, file = paste(projectname, '_qcmat', '.txt', sep=''),
            sep = '\t', quote = FALSE, row.names = FALSE)



##### Running CODEX2


# Estimating library size factor for each sample

Y.nonzero <- Y_qc[apply(Y_qc, 1, function(x){!any(x==0)}),]
pseudo.sample <- apply(Y.nonzero,1,function(x){exp(1/length(x)*sum(log(x)))})
N <- apply(apply(Y.nonzero, 2, function(x){x/pseudo.sample}), 2, median)

save.image(file = paste('CODEX2',"_counts_image.RData", sep=''))




#### Genome-wide normalization using normalize_null ####

#normObj.null <- normalize_null(Y_qc = Y_qc, gc_qc = gc_qc, K = 1:4, N=N)
load('CODEX2_normObj_null.RData')

print("hola")

z.codex <- log(Y_qc/normObj.null$Yhat[[which.max(normObj.null$BIC)]])
cnv_index1 <- which(apply(z.codex,1,sd)>=0.25)


source('/home/proyectos/bioinfo/lodela/VariantCallingFJD/utilities/normalize_codex2_nrl.R') 
normObj_r <- normalize_codex2_nrl(Y_qc = Y_qc, gc_qc = gc_qc, cnv_index = cnv_index1, K = 1:4, N = N)

save(normObj_r, file = paste('CODEX2',"_normObj_r.RData", sep=''))


## Genome-wide normalization with negative control samples

normObj_s <- normalize_codex2_ns(Y_qc = Y_qc,
                               gc_qc = gc_qc, 
                               K = 1:4, 
                               norm_index = 1,
                               N = N)

save(normObj_s , file = paste('CODEX2',"_normObj_s.RData", sep=''))

