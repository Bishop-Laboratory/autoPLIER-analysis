library(PLIER)
data(bloodCellMarkersIRISDMAP)
data(canonicalPathways)
library(tidyverse)
library(rhdf5)

allPaths <- combinePaths(bloodCellMarkersIRISDMAP, canonicalPathways)


h5file <- "data/CellO_data/bulk_RNA_seq_training_set/bulk_log_tpm.h5"
mat <- h5read(h5file, name = "expression")
colnames <- h5read(h5file, name = "experiment")
colnames(mat) <- colnames
genenames <- h5read(h5file, name = "gene_id")
rownames(mat) <- genenames  

library(EnsDb.Hsapiens.v86)
res <- AnnotationDbi::select(x = EnsDb.Hsapiens.v86, keys = as.character(genenames), 
                             columns = c("SYMBOL"))
res <- distinct(res, SYMBOL, .keep_all = TRUE) 

mat2 <- mat[res$GENEID,]
rownames(mat2) <- res$SYMBOL
mat3 <- mat2
rmrow <- rowSums(mat3)
rmrow <- which(rmrow == 0)
mat4 <- mat3[-rmrow,]
rm(mat, mat2, mat3)
gc()
plierResult <- PLIER(mat4, canonicalPathways, seed = 42)
save(plierResult, file = "analysis/PLIER_initial_testing/plierResult.rda")


