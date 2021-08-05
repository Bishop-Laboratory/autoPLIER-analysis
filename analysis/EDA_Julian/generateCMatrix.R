#' @title Generate C Matrix
#' 
#' @description Generates the C Matrix for PLIER
#' 
#' @param category Which categories from msigdbr to include? Should be a 
#' list of vectors containing msigbr categories.
#' @param gene_type Can be "gene_symbol", "ensembl_gene", "entrez_gene"
#' 
#' @importFrom rlang %>% .data !!
#' 
#' @examples 
#' 
#' generateCMatrix(category=list("H" = NA, "C2" = NA, "C5" = c("GO:BP", "GO:MF")), gene_type="ensembl_gene") 
#' 
#' generateCMatrix(category=list("H" = NA, "C2" = NA, "C5" = c("GO:BP", "GO:MF")), gene_type="gene_symbol") 
#' 
#' generateCMatrix(category=list("H" = NA, "C2" = NA, "C5" = c("GO:BP", "GO:MF")), gene_type="entrez_gene") 
#' 
generateCMatrix <- function(category=list("H" = NA,
                                          "C2" = NA, 
                                          "C5" = c("GO:BP")),
                            gene_type="ensembl_gene") {
  
  lapply(names(category), function(cat) {
    subcat <- category[[cat]]
    
    if (! is.na(subcat[1])) {
      lapply(subcat, function(sub_cat_now) {
        msigdbr::msigdbr(category = cat, subcategory = sub_cat_now,
                         species = "Homo sapiens")
      }) %>%
        dplyr::bind_rows()
    } else {
      msigdbr::msigdbr(category = cat, 
                       species = "Homo sapiens")
    }
  }) %>% 
    dplyr::bind_rows() %>%
    dplyr::mutate(yesno=1) %>%
    dplyr::select(gs_name, !! gene_type, yesno) %>%
    dplyr::distinct() %>%
    tidyr::pivot_wider(id_cols = !! gene_type,
                       names_from = gs_name,
                       values_from = yesno,
                       values_fill = 0) %>%
    tibble::column_to_rownames(var = gene_type) %>%
    as.matrix() %>%
    return()
}
