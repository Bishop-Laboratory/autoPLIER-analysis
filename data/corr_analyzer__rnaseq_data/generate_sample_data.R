library(ontologyIndex)
library(GEOmetadb)
library(rhdf5)
library(jsonlite)
setwd(file.path(rprojroot::find_rstudio_root_file(), "data/corr_analyzer__rnaseq_data/"))
source("helper.R")

# Get REGEX dictionary
tissueDict <- jsonlite::read_json("tissueDictionary.json")

# Get the data
dir.create("data/ARCHS4_Download", recursive = TRUE, showWarnings = FALSE)
if (! file.exists("data/ARCHS4_Download/human_tpm_v11.h5")){
  # Download TPM files
  humanTxURL <- "https://s3.amazonaws.com/mssm-seq-matrix/human_tpm_v11.h5"
  system("aws s3 cp --no-sign-request s3://mssm-seq-matrix/human_tpm_v11.h5 data/ARCHS4_Download/human_tpm_v11.h5")
}

if (! file.exists("data/cellosaurus__v40.txt")) {
  download.file("https://ftp.expasy.org/databases/cellosaurus/cellosaurus.txt", "data/cellosaurus__v40.txt")
}

dataFile <- "data/ARCHS4_Download/human_tpm_v11.h5"


if (! file.exists("data/unfiltered_sample_data.rda")) {
  # Get metadata
  sample_id <- h5read(dataFile, "meta/samples/geo_accession")
  tissue <- h5read(dataFile, name = "meta/samples/source_name_ch1")
  series <- h5read(dataFile, name = "meta/samples/series_id")
  organism <- h5read(dataFile, name = "meta/samples/organism_ch1")
  institute <- h5read(dataFile, name = "meta/samples/contact_institute")
  char <- h5read(dataFile, "meta/samples/characteristics_ch1")
  scprob <- h5read(dataFile, "meta/samples/singlecellprobability")
  reads_aligned <- h5read(dataFile, name = "meta/samples/readsaligned")
  Sample_submission_date <- h5read(dataFile, name = "meta/samples/last_update_date")
  reads_total <- h5read(dataFile, name = "meta/samples/readstotal")
  
  sample_data <- data.frame(sample_id, tissue, series, institute, char,
                            Sample_submission_date, 
                            reads_aligned, reads_total)
  
  if (! file.exists("GEOmetadb.sqlite") ) {
    sqlfile <- getSQLiteFile()
  } else {
    sqlfile <- "GEOmetadb.sqlite"
  }
  con <- dbConnect(SQLite(), "GEOmetadb.sqlite")
  metaTable <- dbGetQuery(con,'select gsm,extract_protocol_ch1,extract_protocol_ch2,description,data_processing,title from gsm')
  dbDisconnect(con)
  metaTableNow <- metaTable[metaTable$gsm %in% sample_data$sample_id,]
  sample_data <- merge(x = sample_data, y = metaTableNow, by.x = "sample_id", by.y = "gsm")
  sample_data$percent_aln <- sample_data$reads_aligned/sample_data$reads_total * 100
  sample_data$sample_name <- paste0(sample_data$sample_id, ": ", sample_data$title)
  sample_data$type <- "Unknown"
  sample_data$type[grep(sample_data$char, pattern = "cell line|cell type", perl = TRUE)] <- "cell line"
  unfiltered_sample_data <- sample_data
  save(unfiltered_sample_data, file = "data/unfiltered_sample_data.rda")
} else {
  load("data/unfiltered_sample_data.rda")
}

if (! file.exists("data/filtered_sample_data.rda")) {
  ## Filter sample_data ##
  sample_data <- unfiltered_sample_data
  
  # Reads aligned > 5 million < 100 million
  sample_data <- sample_data[sample_data$reads_aligned > 5e6,]
  sample_data <- sample_data[sample_data$reads_aligned < 1e8,]
  
  # pct reads aln > 60% and <= 96%
  sample_data <- sample_data[sample_data$percent_aln <= 96,]
  sample_data <- sample_data[sample_data$percent_aln > 60,]
  
  # Remove single-cell
  patternSingle <- unlist(tissueDict$SingleCell)
  sample_data <- sample_data[unique(grep(pattern = paste(patternSingle, collapse="|"), invert = T,
                                         x = sample_data$extract_protocol_ch1, ignore.case = T, perl = TRUE)),]
  sample_data <- sample_data[unique(grep(pattern = paste(patternSingle, collapse="|"), invert = T,
                                         x = sample_data$extract_protocol_ch2, ignore.case = T, perl = TRUE)),]
  sample_data <- sample_data[unique(grep(pattern = paste(patternSingle, collapse="|"), invert = T,
                                         x = sample_data$data_processing, ignore.case = T, perl = TRUE)),]
  sample_data <- sample_data[unique(grep(pattern = paste(patternSingle, collapse="|"), invert = T,
                                         x = sample_data$tissue, ignore.case = T, perl = TRUE)),]
  sample_data <- sample_data[unique(grep(pattern = paste(patternSingle, collapse="|"), invert = T,
                                         x = sample_data$description, ignore.case = T, perl = TRUE)),]
  
  # Manually assign disease and tissue types
  sample_data$disease <- "Unknown"
  
  # Remove known problematic studies
  toRemove <- c("GSE46665", "GSE45326", "GSE50957", "GSE68086", "GSE72420Xx-xXGSE72509",
                "GSE88850", "GSE100127", "GSE83402", "GSE68799", "GSE47774Xx-xXGSE47792")
  sample_data <- sample_data[! sample_data$series %in% toRemove,]
  
  # Remove 4su seq
  sample_data <- sample_data[unique(grep(pattern = "\\b4su", invert = T,
                                         x = sample_data$tissue, ignore.case = T, perl = TRUE)),]
  
  # Remove human refrence
  sample_data <- sample_data[unique(grep(pattern = "\\buhrr\\b|\\bbhrr\\b|\\bhbrr\\b|\\burr\\b|Universal", invert = T,
                                         x = sample_data$tissue, ignore.case = T, perl = TRUE)),]
  
  # Remove other unwanted types
  sample_data <- sample_data[unique(grep(pattern = "thrombocyt|extracellular|platelet|exosome", invert = T,
                                         x = sample_data$tissue, ignore.case = T, perl = TRUE)),]
  
  sample_data_filtered <- sample_data
  save(sample_data_filtered, file = "data/filtered_sample_data.rda")
} else {
  load("data/filtered_sample_data.rda")
}


### Wrangle Cellosaurus ###
if (! file.exists("data/cellosaurus_v40_attribute_df_and_id_syn_map.rda")) {
  lines <- readLines("data/cellosaurus__v40.txt", encoding = "UTF-8")
  lines <- paste0(lines, rep("__", length(lines)))
  lines <- paste0(lines, collapse = "")
  res <- stringi::stri_split_fixed(lines, pattern = "__//__")
  res_list <- pbapply::pblapply(seq(length(res[[1]])), function(i){
    entry_now <- res[[1]][i]
    if (entry_now == "") {
      return(NULL)
    }
    sptl <- stringi::stri_split_fixed(entry_now, pattern = "__")
    split_df <- as.data.frame(stringi::stri_split_fixed(sptl[[1]], pattern = "   ", simplify = TRUE))
    split_df$ID <- split_df$V2[1]
    colnames(split_df)[1:2] <- c("attribute", "value")
    split_df
  })
  cello_df <- data.table::rbindlist(res_list, fill = TRUE)
  cello_df <- cello_df[cello_df$ID != "" ,c("ID", "attribute", "value")]
  cello_names_df <- cello_df[cello_df$attribute == "SY",]
  cello_names_df$value <- paste0(cello_names_df$value, "; ", cello_names_df$ID)
  cello_names_delim <- cello_names_df$value
  names(cello_names_delim) <- cello_names_df$ID
  cello_names_undelim <- stringi::stri_split_fixed(cello_names_delim, pattern = "; ", simplify = FALSE)
  cello_names_undelim <- lapply(cello_names_undelim, as.data.frame, drop = FALSE)
  names(cello_names_undelim) <- names(cello_names_delim)
  cello_names_undelim_df <- data.table::rbindlist(cello_names_undelim, use.names = TRUE, idcol = "Name")
  cello_id_syn_map <- cello_names_undelim_df
  colnames(cello_id_syn_map) <- c("ID", "SY")
  save(cello_id_syn_map, cello_df, file = "data/cellosaurus_v40_attribute_df_and_id_syn_map.rda")
} else {
  load("data/cellosaurus_v40_attribute_df_and_id_syn_map.rda")
}

# Only human entries kept
cello_ids_hs <- cello_df$ID[which(cello_df$value == "NCBI_TaxID=9606; ! Homo sapiens")]
cello_df_hs <- cello_df[cello_df$ID %in% cello_ids_hs,]
cello_id_syn_map_hs <- cello_id_syn_map[cello_id_syn_map$ID %in% cello_ids_hs,]
###########################

## categorize samples ##
if (! file.exists("data/filtered_sample_data_init_class.rda")) {
  sample_data <- sample_data_filtered
  sample_data$class_confidence <- NA
  dir.create("data/classification", showWarnings = FALSE)
  
  ### Cell line/type provided ###
  gcl <- sample_data$sample_id[grep(sample_data$char, pattern = "cell line:", perl = TRUE)]
  gct <- sample_data$sample_id[grep(sample_data$char, pattern = "cell type:", perl = TRUE)]
  sample_data_cells <- sample_data[sample_data$sample_id %in% unique(c(gct, gcl)),]
  # Cell line
  sample_data_cells$cell_line <- "Unknown"
  cell_lines <- gsub(sample_data_cells$char, pattern = ".*cell line: (.+)", replacement = "\\1", perl = TRUE)
  cell_lines <- gsub(cell_lines, pattern = "Xx-xX.+", replacement = "", perl = TRUE)
  names(cell_lines) <- sample_data_cells$sample_id
  sample_data_cells$cell_line_raw <- "Unknown"
  sample_data_cells$cell_line_raw[sample_data_cells$sample_id %in% gcl] <- cell_lines[names(cell_lines) %in% gcl]
  # Cell type
  sample_data_cells$cell_type <- "Unknown"
  cell_types <- gsub(sample_data_cells$char, pattern = ".*cell type: (.+)", replacement = "\\1", perl = TRUE)
  cell_types <- gsub(cell_types, pattern = "Xx-xX.+", replacement = "", perl = TRUE)
  names(cell_types) <- sample_data_cells$sample_id
  sample_data_cells$cell_type_raw <- "Unknown"
  sample_data_cells$cell_type_raw[sample_data_cells$sample_id %in% gct] <- cell_types[names(cell_types) %in% gct]
  
  ## Match cell line info directly to cellosaurus ##
  
  # Convert all to lower case
  cello_id_syn_map_hs$SY <- tolower(cello_id_syn_map_hs$SY)
  cello_id_syn_map_hs <- unique(cello_id_syn_map_hs)
  sample_data_cells$cell_line_raw <- tolower(sample_data_cells$cell_line_raw)
  sample_data_cells$cell_type_raw <- tolower(sample_data_cells$cell_type_raw)
  
  # Match strings & merge
  cell_line_matches <- sample_data_cells[sample_data_cells$cell_line_raw %in% cello_id_syn_map_hs$SY,]
  cell_line_matches <- merge(x = cell_line_matches, y = cello_id_syn_map_hs, by.x = "cell_line_raw", by.y = "SY")
  cell_type_matches <- sample_data_cells[sample_data_cells$cell_type_raw %in% cello_id_syn_map_hs$SY,]
  cell_type_matches <- merge(x = cell_type_matches, y = cello_id_syn_map_hs, by.x = "cell_type_raw", by.y = "SY")
  cell_type_matches <- cell_type_matches[,order(match(colnames(cell_type_matches), colnames(cell_line_matches)))]
  cells_matched <- rbind(cell_line_matches, cell_type_matches)
  
  # TODO: Duplicates were not specifically assigned. These need to be handeled later
  cells_matched_dups <- cells_matched[cells_matched$sample_id %in% cells_matched$sample_id[duplicated(cells_matched$sample_id)],]
  
  # CASE: sample from cell line. Perfect match with cellosaurus -- highest confidence level = 5
  cells_matched <- cells_matched[! cells_matched$sample_id %in% cells_matched$sample_id[duplicated(cells_matched$sample_id)],]
  cells_matched$class_confidence <- 5
  save(cells_matched, file = "data/classification/cells_matched.rda")
  
  ## Match cell line info to cellosaurus -- no special characters ##
  sample_data_cells <- sample_data_cells[! sample_data_cells$sample_id %in% cells_matched$sample_id &
                                           ! sample_data_cells$sample_id %in% cells_matched_dups$sample_id,]
  cello_id_syn_map_hs$SY_plain <- tolower(gsub(cello_id_syn_map_hs$SY, pattern = "#|\\[|\\]|\\(|\\)|-|_| |/|\\+", replacement = "", perl = TRUE))
  cello_id_syn_map_hs_alt <- unique(cello_id_syn_map_hs[,c(1,3)])
  cell_line_alt_matches <- sample_data_cells[sample_data_cells$cell_line_raw %in% cello_id_syn_map_hs_alt$SY_plain,]
  cell_line_alt_matches <- merge(x = cell_line_alt_matches, y = cello_id_syn_map_hs_alt, by.x = "cell_line_raw", by.y = "SY_plain")
  cell_type_alt_matches <- sample_data_cells[sample_data_cells$cell_type_raw %in% cello_id_syn_map_hs_alt$SY_plain,]
  cell_type_alt_matches <- merge(x = cell_type_alt_matches, y = cello_id_syn_map_hs_alt, by.x = "cell_type_raw", by.y = "SY_plain")
  cell_type_alt_matches <- cell_type_alt_matches[,order(match(colnames(cell_type_alt_matches), colnames(cell_line_alt_matches)))]
  cells_alt_matched <- rbind(cell_line_alt_matches, cell_type_alt_matches)
  
  # TODO: Duplicates were not specifically assigned. These need to be handeled later
  cells_alt_matched_dups <- cells_alt_matched[cells_alt_matched$sample_id %in% cells_alt_matched$sample_id[duplicated(cells_alt_matched$sample_id)],]
  
  # CASE: sample from cell line. Perfect match with cellosaurus -- highest confidence level = 5
  cells_alt_matched <- cells_alt_matched[! cells_alt_matched$sample_id %in% cells_alt_matched$sample_id[duplicated(cells_alt_matched$sample_id)],]
  cells_alt_matched$class_confidence <- 4 # TODO: Issues with th1 being assigned to CVCL_8306
  save(cells_alt_matched, file = "data/classification/cells_alt_matched.rda")
  
  
  # 
  # #### MetaSRA classification ####
  # if (! file.exists("data/classification/meta_sra_v1-6.json")) {
  #   download.file("http://metasra.biostat.wisc.edu/static/metasra_versions/v1.6/metasra.v1-6.json", 
  #                 destfile = "data/classification/meta_sra_v1-6.json")
  # }
  # ###############################
  # 
  # #### GEOracle classification ####
  # system("cd data/classification && git clone https://github.com/VCCRI/GEOracle.git")
  # #################################
  # 
  # #### URSA-HD classification ####
  # "https://www.sciencedirect.com/science/article/pii/S240547121830509X"
  # ################################
  # 
  # #### CREEDS classification ####
  # 
  # ###############################
  # 
  # #### xCell classification ####
  # # Does ssGSEA for typing bulk tissues -- needs TPM 
  # library(xCell)
  # 
  # ##############################
  # 
  # #### ALE classification ####
  # # Need to run in python -- also need to pip install wrenlab, metalearn, and click
  # system("cd data/classification && git clone https://github.com/wrenlab/label-extraction.git")
  # system("pip install wrenlab metalearn click")
  # ############################
  # 
  # #### phenopredict ####
  # # Very promising -- has integration with recount2
  # # https://pubmed.ncbi.nlm.nih.gov/29514223/
  # 
  # #####################
  # 
  # #### GREIN ####
  # # Not very helpful, but code showing how you can set up
  # # ontology search with metaSRA: https://github.com/uc-bd2k/GREIN/blob/master/R/after_processing.R
  # ###############
  # 
  # #### BioPortal ####
  # # Using annotator plus -- we can do an ontology search from GEO metadata
  # # Pretty useful option particularly for typing samples by disease -- also they have a REST API
  # 
  # ###################
  # 
  
  
  
  # CASE: no direct match made, but successfully matched without special characters
  
  
  
  # Grep on cell line + type using cellosaurus
  sample_data_cells$cell_str <- paste(sample_data_cells$cell_line_raw, sample_data_cells$cell_type_raw)  
  
  
  
  # 
  # # Get list of all cancer cell lines
  # cell_name <- c()
  # cell_alias <- c()
  # bto_vec <- c()
  # disease_vec <- c()
  # 
  # for (i in 1:length(lines)) {
  #   line <- lines[i]
  #   
  #   if (line == "//") {
  #     
  #   }
  # }
  
  
  ## Neither provided ##
  
  
  
  # Prep for grep
  sample_data$char <- gsub(sample_data$char, pattern = "Xx-xX", replacement = " ", perl = T, ignore.case = F)
  newTissue <- paste0(sample_data$tissue, " ", sample_data$char)
  
  # Spot fix
  sample_data$tissue[sample_data$series == "GSE48865"] <- "glioma"
  sample_data$tissue[sample_data$series == "GSE72790"] <- "leukemia"
  sample_data$tissue[sample_data$series %in% c("GSE49155", "GSE39121", "GSE62118")] <- "lung cancer"
  sample_data$tissue[sample_data$series %in% c("GSE65185Xx-xXGSE65186", "GSE78220")] <- "melanoma"
  sample_data$tissue[sample_data$series %in% c("GSE79492", "GSE56066")] <- "breast cancer"
  sample_data$tissue[sample_data$series == "GSE60052" & sample_data$tissue == "tumor"] <- "lung cancer"
  sample_data$tissue[sample_data$series %in% c("GSE47944Xx-xXGSE47965")] <- "skin"
  sample_data$tissue[sample_data$series == "GSE66301"] <- "breast cancer"
  sample_data$tissue[sample_data$series == "GSE66306"] <- "pbmc"
  
  # Re-format tissue names so we can grep cell types
  newTissue = gsub(pattern = "\\.", replacement = " ", x = newTissue)
  newTissue = gsub(pattern = "'", replacement = " ", x = newTissue)
  newTissue = gsub(pattern = "\\(", replacement = " ", x = newTissue)
  newTissue = gsub(pattern = "\\)", replacement = " ", x = newTissue)
  newTissue = gsub(pattern = "\\[", replacement = " ", x = newTissue)
  newTissue = gsub(pattern = "\\]", replacement = " ", x = newTissue)
  newTissue = tolower(newTissue)
  sample_data$tissue <- newTissue
  
  # Categorize cancer samples
  cancerString <- unlist(tissueDict$tumors[[1]])
  normalString <- unlist(tissueDict$tumors[[2]])
  cancerGrep <- grep(sample_data$tissue, pattern = paste0(cancerString, 
                                                          collapse = "|"), 
                     perl = TRUE, ignore.case = TRUE)
  sample_data$disease[cancerGrep] <- "Cancer"
  normalGrep <- grep(sample_data$tissue, pattern = paste0(normalString, 
                                                          collapse = "|"), 
                     perl = TRUE, ignore.case = TRUE)
  sample_data$disease[normalGrep] <- "Unknown"
  
  # hidden tumor samps -- poorly annotated...
  hiddenSamps <- paste0("GSM", c(
    1185643:1185648
  ))
  sample_data$disease[sample_data$sample_id %in% hiddenSamps] <- "Cancer"
  
  # Categorize samples by tissue type using REGEX dictionary
  tissueDictNow <- tissueDict[which(! names(tissueDict) %in% c("tumors", "SingleCell", "stem-like"))]
  sample_data_cat <- categorizeMetaData(metadata = sample_data,
                                        cols = "tissue",
                                        dictionary = tissueDictNow)
  n <- length(colnames(sample_data)) + 1
  m <- n + length(tissueDictNow) - 1 
  rSums <- rowSums(sample_data_cat[,c(n:m)])
  sample_data_cat_main <- unique(sample_data_cat[rSums == 1,]) # Keep only unique assigments
  sample_data_cat_main$category <- "category"
  possibles <- colnames(sample_data_cat_main[,c(n:m)])
  resVec <- pbapply::pbsapply(seq(length(sample_data_cat_main$sample_id)), function(i) {
    possibles[which(sample_data_cat_main[i,c(n:m)] == 1)]
  })
  sample_data_cat_main$category <- resVec
  save(sample_data_cat, rSums, n, m, sample_data, tissueDictNow, sample_data_cat_main, file = "data/filtered_sample_data_init_class.rda")
} else {
  load("data/filtered_sample_data_init_class.rda")
}


# Samples which received no classification -- use Cellosaurus for secondary classification
sample_data_cat_no_class <- unique(sample_data_cat_main[rSums == 0,]) 
# -- Get list of cancer cells and cancer types from cellosaurus
if (! file.exists("data/cellDisease.rda")) {
  lines <- readLines("data/cellosaurus__v40.txt", encoding = "UTF-8")
  
  # lines <- paste0(lines, collapse = "")
  
  # Get list of all cancer cell lines
  cell_name <- c()
  cell_alias <- c()
  bto_vec <- c()
  diseaseVec <- c()
  cancerCellVec <- c()
  
  for (i in 1:length(lines)) {
    if (i %% 1000 == 0) { message(i) }
    line <- lines[i]
    grepRes1 <- grep(line, pattern = "^ID")
    grepRes2 <- grep(line, pattern = "^SY")
    grepRes3 <- grep(line, pattern = "^CA")
    grepRes4 <- grep(line, pattern = "^DI")
    if (length(grepRes1)) {
      tumors <- FALSE
      disease <- NULL
      cellsNow <- gsub(x = line, pattern = "^ID\\s+",
                       replacement = "", perl = TRUE)
    }
    if (length(grepRes2)) {
      cellsNow2 <- gsub(x = line, pattern = "^SY\\s+",
                        replacement = "", perl = TRUE)
      cellsNow <- c(cellsNow, unlist(strsplit(cellsNow2, split = "; ")))
    }
    if (length(grepRes4)) {
      disease <- gsub(x = line, pattern = "^DI.+;.+; ",
                      replacement = "", perl = TRUE)
    }
    if (length(grepRes3)) {
      if (line == "CA   Cancer cell line") {
        if (is.null(disease)) {
          next
        }
        cancerCellVec <- c(cancerCellVec, cellsNow)
        diseaseVec <- c(diseaseVec, rep(disease, length(c(cellsNow))))
        if (! length(cancerCellVec) == length(diseaseVec)) {
          stop()
        }
      }
    }
  }
  
  cellDisease <- data.frame(
    "cell" = cancerCellVec,
    "disease" = diseaseVec
  )
  save(cellDisease, file = "data/cellDisease.rda")
  
} else {
  load("data/cellDisease.rda")
}
# -- grep cell types and assign disease
if (! file.exists("data/tissueDisease.rda")) {
  # Clean for grep
  cellDisease$cell <- gsub(cellDisease$cell, pattern = "\\+", replacement = "\\\\+")
  cellDisease$cell <- gsub(cellDisease$cell, pattern = "\\[", replacement = " ")
  cellDisease$cell <- gsub(cellDisease$cell, pattern = "\\]", replacement = " ")
  cellDisease$cell <- gsub(cellDisease$cell, pattern = "\\(", replacement = " ")
  cellDisease$cell <- gsub(cellDisease$cell, pattern = "\\)", replacement = " ")
  cellDisease$cell <- gsub(cellDisease$cell, pattern = "\\.", replacement = "\\\\.")
  cellDisease$cell <- gsub(cellDisease$cell, pattern = "\\^", replacement = "\\\\^")
  cellDisease$cell <- gsub(cellDisease$cell, pattern = "\\$", replacement = "\\\\$")
  cellDisease$cell <- gsub(cellDisease$cell, pattern = "\\?", replacement = "\\\\?")
  cellDisease$cell <- gsub(cellDisease$cell, pattern = "\\|", replacement = "\\\\|")
  cellDisease$cell <- gsub(cellDisease$cell, pattern = "\\*", replacement = "\\\\*")
  cellDisease$cell2 <- gsub(cellDisease$cell, pattern = "\\#", replacement = "")
  
  # Prepare for vec
  tissueUnique <- unique(sample_data_cat_no_class$tissue)
  
  # Pattern to match cell types with or without add info separated by - or _
  # This basically says, there can be something such that "-cellName-", 
  # "_cellName_" or " cellName " and all the possible combinations.
  cellDisease$cell1 <- sapply(cellDisease$cell, FUN = function(x) {
    paste0("\\b", x, "[_-]+|\\b", x, "\\b|[_-]+", x, "\\b|[_-]+",x,"[_-]+")
  })
  cellDisease$nchar <- sapply(cellDisease$cell2, FUN = nchar)
  cellDisease <- cellDisease[cellDisease$nchar > 2,]
  
  tableDF <- as.data.frame(table(cellDisease$disease))
  cellDiseaseList <- split(cellDisease, f = cellDisease$disease)
  calculateGrep <- function(wordVec, patternVec) {
    # wordVec <-  tissueUnique[1:100]
    # patternVec <- cancerCellVec1[30344:3344]
    return(sapply(wordVec, FUN = grep, ignore.case = TRUE,
                  pattern = paste0(as.character(patternVec), collapse = "|"), perl = TRUE))
  }
  
  # Split the tissueUnique -- from https://stackoverflow.com/questions/3318333/split-a-vector-into-chunks-in-r
  chunk2 <- function(x,n) split(x, cut(seq_along(x), n, labels = FALSE)) 
  tissueUniqueList <- chunk2(tissueUnique, 100)
  
  # apply to each chunk -- map GEO tissue annotations to cellosaurus diseases
  tissueDiseaseList <- parallel::mclapply(
    tissueUniqueList, 
    mc.cores = length(tissueUniqueList),
    FUN = function(tissueNow) {
      namesList <- c()
      diseaseList <- c()
      for (i in 1:length(cellDiseaseList)) {
        print(i)
        dataNow <- cellDiseaseList[[i]]
        diseaseNow <- names(cellDiseaseList)[i]
        cancerVecNow <- unique(dataNow$cell1)
        cancerVecNowSplit <- split(cancerVecNow, ceiling(seq_along(cancerVecNow)/100)) # Chunks of 500
        names(cancerVecNowSplit) <- paste0(names(cancerVecNowSplit), "_xxx_")
        grepRes <- lapply(cancerVecNowSplit, FUN = function(splitNow) {
          tryCatch({calculateGrep(tissueNow, patternVec = splitNow)},
                   error = function(e) {})
        })
        grepRes <- unlist(grepRes, recursive = T, use.names = T)
        namesNow <- gsub(names(grepRes), pattern = ".+_xxx_\\.", replacement = "")
        namesList <- c(namesList, namesNow)
        diseaseList <- c(diseaseList, rep(diseaseNow, length(namesNow)))
        if (! length(namesList) == length(diseaseList)) {
          stop()
        }
      }
      
      tissueDisease <- data.frame(
        "tissue" = namesList,
        "disease" = diseaseList
      )
      tissueDisease
    }
  )
  tissueDisease <- data.table::rbindlist(tissueDiseaseList)
  ## Extra cleanup -- some inaccurate cellosaurus entries... ##
  disease2 <- sapply(tissueDisease$disease, FUN = function(x) {
    paste0("\\b", x, "[_-]+|\\b", x, "\\b|[_-]+", x, "\\b|[_-]+",x,"[_-]+")
  })
  names(disease2) <- tissueDisease$disease
  disease2 <- disease2[! duplicated(disease2)]
  tissueDisease$alt_tissue <- "None"
  for (i in 1:length(disease2)) {
    hits <- grep(tissueDisease$tissue, pattern = disease2[i], perl = TRUE, ignore.case = TRUE)
    tissueDisease$alt_tissue[hits][which(tissueDisease$disease[hits] != names(disease2)[i])] <- names(disease2)[i]
  }
  tissueDisease$disease[tissueDisease$alt_tissue != "None"] <- tissueDisease$alt_tissue[tissueDisease$alt_tissue != "None"]
  tissueDisease <- tissueDisease[,c(-3)]
  
  # Remove non-human classes
  tissueDisease <- tissueDisease[grep(tissueDisease$disease, pattern = "\\bhamster\\b|\\brat\\b|\\bmouse\\b|\\bcanine\\b|\\bfeline\\b|Li-Fraumeni|Mycosis|\\bgoldfish\\b",
                                      ignore.case = T, invert = T),]
  
  # Clean up duplicates
  tissueDisease <- tissueDisease[which(! duplicated(tissueDisease$tissue)),]
  
  save(tissueDisease, file = "data/tissueDisease.rda")
} else {
  load("data/tissueDisease.rda")
}


if (! file.exists("data/sample_data_raw.rda")) {
  # Wrangle
  sample_data_cat_cell <- sample_data_cat_no_class[sample_data_cat_no_class$tissue %in% tissueDisease$tissue,]
  sample_data_cat_cell <- merge(x = sample_data_cat_cell[,c(1:(n-3))], y = tissueDisease, by = "tissue")
  sample_data_cat_cell$tissue_orig <- sample_data_cat_cell$tissue
  sample_data_cat_cell$tissue <- sample_data_cat_cell$disease
  sample_data_cat_cell$disease <- "Cancer"
  
  # Make final tissue assignments
  sample_data_cat_cell_cat <- categorizeMetaData(metadata = sample_data_cat_cell,
                                                 cols = "tissue",
                                                 dictionary = tissueDictNow)
  
  n <- length(colnames(sample_data))
  m <- n + length(tissueDictNow) - 2
  
  rSums <- rowSums(sample_data_cat_cell_cat[,c((n+1):(m+1))])
  sample_data_cat_cell_cat_main <- unique(sample_data_cat_cell_cat[rSums == 1,]) # Keep only unique assigments
  sample_data_cat_cell_cat_main$category <- "category"
  possibles <- colnames(sample_data_cat_cell_cat_main[,c((n+1):(m+1))])
  resVec <- pbapply::pbsapply(seq(sample_data_cat_cell_cat_main$sample_id), function(i) {
    possibles[which(sample_data_cat_cell_cat_main[i,c((n+1):(m+1))] == 1)]
  })
  sample_data_cat_cell_cat_main$category <- resVec
  
  # Combine data
  sample_data_orig <- sample_data
  sample_data_cat_cell_cat_main <- sample_data_cat_cell_cat_main[, colnames(sample_data_cat_cell_cat_main) %in% colnames(sample_data_cat_main)]
  sample_data_cat_main <- sample_data_cat_main[,
                                               colnames(sample_data_cat_main) %in% colnames(sample_data_cat_cell_cat_main)
  ]
  all(colnames(sample_data_cat_cell_cat_main) %in% colnames(sample_data_cat_main))
  sample_data_cat_cell_cat_main <- sample_data_cat_cell_cat_main[,
                                                                 order(match(colnames(sample_data_cat_cell_cat_main), colnames(sample_data_cat_main)))
  ]
  all(colnames(sample_data_cat_cell_cat_main) == colnames(sample_data_cat_main))
  
  sample_data_cat <- rbind(sample_data_cat_main, sample_data_cat_cell_cat_main)
  sample_data$category <- "Unknown"
  sample_data_cat <- sample_data_cat[,colnames(sample_data_cat) %in% colnames(sample_data)]
  sample_data <- sample_data[, -which(colnames(sample_data) == "class_confidence")]
  all(colnames(sample_data_cat) == colnames(sample_data))
  sample_data_uncat <- sample_data[! sample_data$sample_id %in% sample_data_cat$sample_id,]
  sample_data <- rbind(sample_data_cat, sample_data_uncat)
  sample_data <- sample_data[! duplicated(sample_data$sample_id),]
  
  # Final clean-up
  sample_data$disease[sample_data$disease == "Unknown"] <- "Normal"
  sample_data$category[sample_data$series %in% c("GSE140442")] <- "stem-like"
  sample_data$disease[sample_data$series %in% c("GSE140442")] <- "Normal"
  sample_data <- sample_data[! is.na(sample_data$category),]
  
  # Cell line vs tissue
  gcl <- grepl(sample_data$char, pattern = "cell line")
  
  # This constitutes the colData prior to dim reduction assisted assignment
  sample_data_raw <- sample_data
  save(sample_data_raw, file = "data/sample_data_raw.rda")
} else {
  load("data/sample_data_raw.rda")
}

library(tidyverse)
sample_data <- sample_data_raw %>% 
  filter(
    disease == "Normal", 
    category != "ewing sarcoma"
  )

## Get the expression to the gene lvl
if (! file.exists("data/genelvl_tpm.rds")) {
  sample_id <- h5read(dataFile, "meta/samples/geo_accession")
  genes <- h5read(dataFile, name = "meta/transcripts/gene_symbol")
  
  # Need to aggregate and compress...
  d <- seq(sample_id)
  sampch <- split(d, ceiling(seq_along(d)/200))
  sampch2 <- split(sampch, ceiling(seq_along(sampch)/44))
  
  pbapply::pblapply(
    seq(sampch2),
    function(j) {
      sampchnow <- sampch2[[j]]
      if (! file.exists(paste0("data/aggregate_tpm/agg__", j, ".rds"))) {
        aggres <- parallel::mclapply(seq(sampchnow), function(i){
          sampsnow <- sampchnow[[i]]
          mat <- h5read(dataFile, name = "data/expression", index = list(sampsnow, NULL))
          toagg <- mat %>%
            t() %>%
            as_tibble()
          colnames(toagg) <- sample_id[sampsnow]
          toagg$gene_name <- genes
          aggfin <- toagg %>% 
            filter(gene_name != "missing") %>%
            group_by(gene_name) %>% 
            summarise(across(.fns = sum))
          if(i != 1) {
            select(aggfin, -gene_name)
          } else {
            aggfin
          }
        }, mc.cores = 44)
        aggfull <- bind_cols(aggres)
        write_rds(aggfull, file = paste0("data/aggregate_tpm/agg__", j, ".rds"))
      }
      
      return(NULL)
    }
  )
  
  reslst <- pbapply::pblapply(
    seq(sampch2),
    function(i) {
      aggfin <- read_rds(file = paste0("data/aggregate_tpm/agg__", i, ".rds"))
      if(i != 1) {
        select(aggfin, -gene_name)
      } else {
        aggfin
      }
    }
  )
  tpmfull <- bind_cols(reslst)
  write_rds(tpmfull, file = "data/genelvl_tpm.rds")
  tpmfull <- readRDS("data/genelvl_tpm.rds")
  tpmsmall <- tpmfull %>% select(gene_name, sample_data$sample_id[sample_data$sample_id %in% colnames(tpmfull)])
  save(tpmsmall, sample_data, file = "data/tpm_normaltissue.rda")
  load("data/tpm_normaltissue.rda")
  sample_datasmall <- sample_data[sample_data$sample_id %in% colnames(tpmsmall),]
  all(colnames(tpmsmall)[-1] == sample_datasmall$sample_id)
  tpmsmall[,1] %>% mutate(gene_name2 = is.list(gene_name)) -> dd
  is_bad_column <- vapply(tpmsmall[,1:10], function(xx) !is.null(dim(xx)), logical(1))
  dim(tpmsmall$gene_name) <- NULL
  write_csv(x = tpmsmall, file = "data/tpm_normaltissue.csv")
  write_csv(x = sample_datasmall, file = "data/tpm_normaltissue__colData.csv")
}



# Downsamples
library(recount)
rse <- recount::download_study(project = "SRP012682", download=FALSE)
# Downloaded data using google chrome through X11
load("data/rse_gene.Rdata")
tpm <- recount::getTPM(rse = rse_gene)

write_tsv(
  as.data.frame(tpm), file = "data/gtex_tpm.tsv"
)
dd <- colData(rse_gene) %>% as.data.frame()

all(rownames(dd) == colnames(tpm))
colData(rse_gene) %>% as.data.frame() %>% write_tsv("data/gtex_tpm_coldata.tsv")






