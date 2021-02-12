# Analysis

This folder will contain the results of various analyses performed in support for the overall project. As new analyses are conceived, they will be added here. 

Each analysis folder will contain a `README` which gives additional details. 

## List of Analyses

The following list will correspond to the analyses which have been planned for this project.

#### Testing PLIER on the CellO training set

Data source: normalized RNA-Seq - `data/CellO_data/bulk_RNA_seq_training_set/bulk_log_tpm.h5`

Folder: `PLIER_initial_testing_Henry/`

Questions:

1. Does PLIER allow for the definition of cell types easily?
2. How much more will we have to do in addition to using PLIER?

#### EDA of CellO training dataset 

Data source: normalized RNA-Seq - `data/CellO_data/bulk_RNA_seq_training_set/bulk_log_tpm.h5`

Justin's folder - `analysis/CellO_training_data_eda_Justin`

Serhiy's folder - `analysis/CellO_training_data_eda_Serhiy`

Recommended Tasks (feel free to change):

1. Load the dataset
2. Scale the dataset if needed (VST transform might be useful to try as well)
3. Generate PCA and choose the number of dimensions to keep (elbow plot helps with this)
4. Do nearest neighbor search and generate neighborhood graph (Seurat/scanpy has builtin functions for this)
5. Do clustering (Seurat/scanpy has built-in functions for this)
6. Do UMAP (Seurat/scanpy has built-in functions for this)
7. Overlay the flat labels (once Matt provides them)

Questions: 
1. How well do the labels match up with the clusters?
2. Does the UMAP visualization represent these labels well?




