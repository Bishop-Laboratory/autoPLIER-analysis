
library(PLIER)
data(canonicalPathways)

print("reading the source matrix")
data_mat <- as.matrix(read.table("mat3.csv", header=T, sep = ","))

seed = 25
num_lvs = 150

print(paste("runing PLIER on new LV dims:", num_lvs))
outputFile = paste0("results/plierResult-train-vm-seed_", seed, "-")
outputFile = paste0(outputFile, num_lvs, "_lvs.rda")
plierResult <- PLIER(data_mat, canonicalPathways, k=num_lvs, seed=seed)
save(plierResult, file = outputFile)


print("DONE")