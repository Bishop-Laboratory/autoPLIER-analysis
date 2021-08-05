# Set R Root
setwd(rprojroot::find_rstudio_root_file())

# Libraries
library(PLIER)
library(tidyverse)
library(rhdf5)
source("analysis/EDA_Julian/generateCMatrix.R")

# MAGIC VALUES
SEED = 25
NUMBER_LVS = 150
TRAIN_DATA_PATH <- "data/CellO_data/bulk_RNA_seq_training_set/bulk_log_tpm.h5"
OUTPUT_PLIER = paste0("results/plierResult-train-vm-seed_", SEED, "-", NUMBER_LVS, "_lvs.rda")
C_MAT <- generateCMatrix(category = list("H" = NA, "C2" = NA, "C5" = c("GO:BP")),
                         gene_type="ensembl_gene")

# Get the data
if (! file.exists(TRAIN_DATA_PATH)) {
  download.file("https://cello-multiplier.s3.us-west-2.amazonaws.com/bulk_log_tpm.h5",
                destfile = TRAIN_DATA_PATH)
}

# Load the data
mat <- h5read(TRAIN_DATA_PATH, name = "expression")
colnames <- h5read(TRAIN_DATA_PATH, name = "experiment")
colnames(mat) <- colnames
genenames <- h5read(TRAIN_DATA_PATH, name = "gene_id")
rownames(mat) <- genenames  

# Wrangle the data
mat2 <- mat
rmrow <- rowSums(mat)
rmrow <- which(rmrow == 0)
mat3 <- mat2[-rmrow,]
print("reading the source matrix")
data_mat <- mat3

# Run PLIER
print(paste("runing PLIER on new LV dims:", NUMBER_LVS))
plierResult <- PLIER(data = data_mat, priorMat = C_MAT, k=NUMBER_LVS, seed=SEED)
save(plierResult, file = OUTPUT_PLIER)


print("DONE")