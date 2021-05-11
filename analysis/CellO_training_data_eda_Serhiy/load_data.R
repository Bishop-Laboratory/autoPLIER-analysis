# loads into plierResult varible
load(file = "data/plierResult-train-v2.rda")

write.table(plierResult$B, file="data/plierResult-cello_train/B.csv", sep = ",")
write.table(plierResult$Z, file="data/plierResult-cello_train/Z.csv", sep = ",")
# plierResult$Z