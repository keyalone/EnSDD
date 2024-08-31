#' This function focuses on cleaning SRT data
#' and store the SRT expression data, spatial location and image path in the Seurat object
#'
#'
#' @param counts_path the .txt path for SRT gene expression data. The file need include gene names and spots names.
#' @param loc_path the .txt path for meta data of the spots.
#' @param img_path the .tif path of H\& image of the corresponding sequencing tissue.
#' @param n_HVG the number of Highly Variation Genes (HVGs) selected by Seurat.
#' @param n_PCA the number of low-dimensional of SRT gene expression.
#'
#' @return a Seurat object.
#'
#'
#' @export
data_process <- function(counts_path, loc_path, img_path, n_HVG = 2000, n_PCA = 20){

  ### input counts file
  counts = read.table(counts_path, sep = '\t')
  counts_ST <- as.matrix(counts)

  ### input location file
  coordinates <- data.frame(read.table(loc_path, sep = '\t'))

  colnames(counts_ST) <- rownames(coordinates)


  ### create Seurat object
  suppressWarnings({Seurat.ST <- Seurat::CreateSeuratObject(counts = counts_ST, assay = "Spatial",
                                                            meta.data = coordinates)})
  Seurat.ST <- Seurat::SCTransform(Seurat.ST, assay = "Spatial", verbose = FALSE)

  ### select HVG genes
  Seurat.ST <- Seurat::FindVariableFeatures(Seurat.ST, nfeatures = n_HVG, verbose = FALSE)
  ### PCA
  Seurat.ST <- Seurat::ScaleData(Seurat.ST, verbose = FALSE)
  Seurat.ST <- Seurat::RunPCA(Seurat.ST, npcs = n_PCA, verbose = FALSE)

  data.list <- list(counts_path = counts_path, loc_path = loc_path, img_path = img_path)
  Seurat.ST@images <- data.list
  return(Seurat.ST)
}

############ transfer clustering results to connectivity matrix ##########
#' This function will transfer the spots label to binary similarity matrix
#'
#' @param  memb the spots label vector
#'
#' @return a binary spots similarity matrix
clust2Mat <- function(memb){
  N <- length(memb)
  mat_cn <- as.numeric(outer(memb, memb, FUN = "==")) - outer(1:N, 1:N, "==")
  rownames(mat_cn) <- names(memb)
  colnames(mat_cn) <- names(memb)
  return(mat_cn)
}


########### ensemble strategy ###############
#' The adaptive weighted ensemble-based learning method to integrate the multiple binary spots similarity matrix
#'
#'
#' @importFrom parallel makeCluster stopCluster parApply
#' @importFrom abind abind
#'
#' @param Results.clustering a list contains all the results of individual similarity matrix. The elements of list is a matrix, spots * spots.
#' @param lambda hyper-parameter constrain the weight of individual methods for ensemble. If the parameter is set to NULL, then, we will adopt the value in our algorithm.
#' @param prob.quantile numeric of probabilities with values in [0,1]. Default setting is 0.5.
#' @param niter a positive integer represents the maximum number of updating algorithm. Default setting is 100.
#' @param epsilon a parameter represents the stop criterion.
#'
#' @return a list contains a matrix of the ensemble similarity of spots and a vector of the weight assigned to base results.
#'
#'@export


solve_ensemble <- function(Results.clustering, lambda = NULL, prob.quantile = 0.5,
                           niter = 100, epsilon = 1e-5){
  options(digits = 7)
  # Results.clustering <- Results.clustering.all[[1]]
  num.methods <- length(Results.clustering)
  num.spots <- nrow(Results.clustering[[1]])
  num.cell.type <- ncol(Results.clustering[[1]])

  ## initialization V by the mean of individual values
  w <- c(rep(1/num.methods, num.methods))
  H <-  Reduce("+", Map("*", Results.clustering, w))

  if(is.null(lambda)){
    cat("We will adpote a value for lambda in our algorithm...", "\n")
  }

  k <- 1

  while (k <= niter) {
    if(k == 1){
      loss_all_temp <- 0
      temp2 <-  sapply(Results.clustering, L2_norm, Y = H)
      if(is.null(lambda)){
        lambda <- quantile(temp2, probs = prob.quantile)
      }
    }else{
      loss_all_temp <- loss_all
    }
    ##### update w
    temp2 <-  sapply(Results.clustering, L2_norm, Y = H)
    w <- exp(-temp2/lambda)/sum(exp(-temp2/lambda))
    H <-  Reduce("+", Map("*", Results.clustering, w))

    loss_main <- sum(sapply(Results.clustering, L2_norm, Y = H) * w)
    loss_entropy <- sum(w * log(w))
    loss_all <- loss_main + lambda * loss_entropy

    # cat("iter: ", k, "loss_main: ", loss_main, "loss_entropy: ", loss_entropy,
    #     "loss_all: ", loss_all, "lambda: ", lambda, "\n")
    if(k == niter)
      cat("The method maybe not convergens, the algorithm need an larger max_epoches!", "\n")

    diff_iter <- abs(loss_all - loss_all_temp)
    cat("iter: ", k, "loss_main: ", loss_main, "loss_entropy: ", loss_entropy,
        "loss_all: ", loss_all, "lambda: ", lambda, "diff",
        diff_iter, "\n")

    if(diff_iter < epsilon | k >= niter)
      break
    k <- k + 1

  }
  colnames(H) <- colnames(Results.clustering[[1]])
  return(list(H = H, w = w))

}


L2_norm <- function(X, Y){
  return(sqrt(sum((X-Y)^2)))
}

weighted_average <- function(x, w){
  index_non_zero <- which(x != 0, arr.ind = TRUE)
  x_use <- x[index_non_zero]
  w_use <- w[index_non_zero]
  results <- sum(x_use * w_use)
  return(results)
}
########### ensemble strategy ###########

########### patch dealing ##############
#' This function aims to integrate the results obtained by user own running
#' spatial domain methods.
#'
#'
#'
#' @param res a list contains all the results of base SDD methods contain in EnSDD methods.
#' @param clustering_df a data frame contains the clustering label generated by SDD methods generated by users. The row represents the spot and the column represents the methods.
#' @return a list contains clustering labels, spot-spot similarity matrix and running time.
patch_deal <- function(res, clustering_df){
  if(is.null(colnames(clustering_df)))
    stop('Please provide the methods of clustering results')
  res.final = res

  n_path <- ncol(clustering_df)
  names_all <- c(names(res.final[[1]]), colnames(clustering_df))
  spot_names = names(res.final[[1]][[1]])
  ori_len = length(res.final[[1]])

  clustering_df = clustering_df[spot_names,]

  ## integrate
  for (i in 1:n_path) {
    tmp = clustering_df[,i]
    names(tmp) = spot_names
    tmp_sim_mat = clust2Mat(memb = tmp)
    res.final[[1]][[ori_len + i]] = tmp
    res.final[[2]][[ori_len + i]] = tmp_sim_mat
  }
  names(res.final[[1]]) = names_all
  names(res.final[[2]]) = names_all
  return(res.final)
}

