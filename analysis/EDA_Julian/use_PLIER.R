library(PLIER)
data(bloodCellMarkersIRISDMAP)
data(canonicalPathways)
library(tidyverse)
library(rhdf5)






h5file <- "bulk_log_tpm.h5"
mat <- h5read(h5file, name = "expression")
colnames <- h5read(h5file, name = "experiment")
colnames(mat) <- colnames
genenames <- h5read(h5file, name = "gene_id")
rownames(mat) <- genenames  

library(EnsDb.Hsapiens.v86)
res <- AnnotationDbi::select(x = EnsDb.Hsapiens.v86, keys = as.character(genenames), 
                             columns = c("SYMBOL"))
res <- distinct(res, SYMBOL, .keep_all = TRUE) 

mat2 <- mat
rmrow <- rowSums(mat)
rmrow <- which(rmrow == 0)
mat3 <- mat2[-rmrow,]
rm(mat, mat2)
gc()

write.table(mat3, file="mat3.csv", sep = ",")
