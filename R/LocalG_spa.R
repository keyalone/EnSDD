#' Local Getis and Ord's Gi for measure the spatial autocorrelation of selected gene
#'
#' @importFrom spdep localG knearneigh knn2nb
#'
#' @param Seurat_data a SeuratObject created by data_process function in EnSDD R package
#' @param gene_name a selected gene for the calculation of Local Getis and Ord's Gi
#' @param k  the number of nearest neighbours to be returned in the KNN
#'
#' @return a vector represents the Local Getis and Ord's Gi value for a selected gene across samples
#'
#' @export


LocalG_spa <- function(Seurat_data, gene_name, k = 6){

  normalized_exp_data <-  as.matrix(Seurat_data@assays$SCT@data)
  df <- Seurat_data@meta.data[, c("x", "y")][colnames(normalized_exp_data),]

  xycoords <- cbind(df$x, df$y)
  require(spdep)
  knn_object <- spdep::knn2nb(spdep::knearneigh(xycoords, k = k))
  w <- spdep::nb2listw(knn_object, style = "W")


  values <- normalized_exp_data[gene_name, ]
  G_value <- spdep::localG(values, w)

  G_value_re_order <- G_value[colnames(normalized_exp_data)]

  return(G_value)

}
