library(PLIER)
data(canonicalPathways)

data_mat <- as.matrix(read.table("data/CellO_data/bulk_RNA_seq_training_set/plier_train.csv", header=T, sep = ","))

plierResult <- PLIER(train_data, canonicalPathways, seed = 42)

save(plierResult, file = "data/plierResult-train.rda")