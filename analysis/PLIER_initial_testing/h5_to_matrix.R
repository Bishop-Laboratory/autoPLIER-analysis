# To test PLIER via the online interface (http://gobie.csb.pitt.edu/PLIER/#)
# Need to convert .h5 log TPM into a count matrix

library(rhdf5)
library(tidyverse)

h5file <- "data/CellO_data/bulk_RNA_seq_training_set/bulk_log_tpm.h5"
h5ls(h5file)  # See the contents

mat <- h5read(h5file, name = "expression")
colnames <- h5read(h5file, name = "experiment")
colnames(mat) <- colnames
genenames <- h5read(h5file, name = "gene_id")
mat <- cbind(genenames, as.data.frame(mat))

write_tsv(mat, file = "analysis/PLIER_initial_testing/bulk_log_tpm.tsv")

