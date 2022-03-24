library(rhdf5)
library(DESeq2)
library(uwot)
library(Seurat)
library(jsonlite)

source("helper.R")


# Download files from ARCHS4
humanExpURL <- "s3://mssm-seq-matrix/human_matrix.h5"
outFile <- "data/ARCHS4_Download/human_matrix.h5"
if (! file.exists(outFile)) {
  system(paste("aws s3 cp", humanExpURL, outFile))
}

expression <- h5read(outFile, "data/expression")
genes <- h5read(outFile, name = "meta/genes")
samples <- h5read(outFile, "meta/samples/geo_accession")


if (! file.exists("Data/vsd_for_corr.rda")) {
  load("Data/fullRawCountsFiltered.rda")
  colDataFinal <- colDataFinal[colDataFinal$samples %in% colnames(expression),]
  colDataFinal <- colDataFinal[order(match(colDataFinal$samples, colnames(expression))),]
  all(colnames(expression) == colDataFinal$samples)
  rownames(colDataFinal) <- colDataFinal$samples
  # Keep genes expressed in 10%+ of samples
  nonZeroCount <- apply(expression, 1, nonZeroSamps)
  keepInd <- which(nonZeroCount > (length(colnames(expression)) * .1))
  expression <- expression[keepInd,]
  
  # Transform and normalize
  dds <- DESeqDataSetFromMatrix(expression, colData = colDataFinal, 
                                design = ~1)
  # from https://support.bioconductor.org/p/62246/#62250
  inds <- rownames(expression)
  geoMeansList <- mclapply(inds, FUN = function(ind) {
    row <- expression[ind,]
    if (all(row == 0)) {
      0
    } else {
      exp(sum(log(row[row != 0]))/length(row))
    }
  }, mc.cores = 25)
  geoMeans <- unlist(geoMeansList)
  dds <- estimateSizeFactors(dds, geoMeans=geoMeans)
  vsd <- vst(dds)
  timestamp()
  cat("\nDone. Saving vst data...\n")
  save(vsd, file = "Data/vsd_for_corr.rda")
} else {
  cat("\nVSD found -- loading...\n")
  load("Data/vsd_for_corr.rda")
}























