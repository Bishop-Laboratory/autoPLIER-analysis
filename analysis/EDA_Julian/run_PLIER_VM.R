
library(PLIER)
data(canonicalPathways)
data(bloodCellMarkersIRISDMAP)
data(canonicalPathways)
library(tidyverse)
library(rhdf5)
library(EnsDb.Hsapiens.v86)


h5file <- "bulk_log_tpm.h5"
mat <- h5read(h5file, name = "expression")
colnames <- h5read(h5file, name = "experiment")
colnames(mat) <- colnames
genenames <- h5read(h5file, name = "gene_id")
rownames(mat) <- genenames  

mat2 <- mat
rmrow <- rowSums(mat)
rmrow <- which(rmrow == 0)
mat3 <- mat2[-rmrow,]
rm(mat, mat2)
gc()

write.table(mat3, file="mat3.csv", sep = ",")


print("reading the source matrix")
data_mat <- mat3

seed = 25
num_lvs = 150

print(paste("runing PLIER on new LV dims:", num_lvs))
outputFile = paste0("results/plierResult-train-vm-seed_", seed, "-")
outputFile = paste0(outputFile, num_lvs, "_lvs.rda")
plierResult <- PLIER(data_mat, canonicalPathways, k=num_lvs, seed=seed)
save(plierResult, file = outputFile)


print("DONE")