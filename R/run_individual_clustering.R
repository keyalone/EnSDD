# .onLoad <- function(libname, pkgname) {
#   user_permission <- utils::askYesNo("Install miniconda? If you have create the EnSDD environment
#                                      using Anaconda and equip with the essential libray, please
#                                      select no,If you don't create the environment,
#                                      please downloads 50MB and takes time")
#
#   if (isTRUE(user_permission)) {
#     reticulate::install_miniconda()
#   } else {
#     message("Enjoy the EnSDD")
#   }
# }

# EnSDD_pyexec <- function(pyfile = NULL) {
#     reticulate::source_python(pyfile)
# }


############# use gene expression and spatial locations #############
# include BayesSpace, DR.SC, SpaGCN, Giotto-HM
#' platform = "Visium" or "ST"
apply_BayesSpace <- function(Seurat.data, platform = "Visium", n.PCs = 15, n.HVGs = 2000,
                             n.setting = 7, nrep = 50000){
  ### extract dataset
  exp_data <- list(counts = as.matrix(Seurat.data@assays$SCT@counts), logcounts = as.matrix(Seurat.data@assays$SCT@data))
  df <- Seurat.data@meta.data[, c("x", "y", "pixel_x", "pixel_y")]
  temp_colnames <- colnames(df)
  temp_colnames <- sub("x", "row", temp_colnames)
  temp_colnames <- sub("y", "col", temp_colnames)
  colnames(df) <- temp_colnames
  sce <- SingleCellExperiment::SingleCellExperiment(assay = exp_data, colData = df)
  dec <- scran::modelGeneVar(sce)
  top <- scran::getTopHVGs(dec, n = n.HVGs)

  set.seed(102)
  sce <- scater::runPCA(sce, subset_row=top)
  sce <- BayesSpace::spatialPreprocess(sce, platform=platform, skip.PCA=TRUE)

  if(platform == "Visium"){
    gamma <- 3
  }
  if(platform == "ST"){
    gamma <- 2
  }

  set.seed(150)
  # if(is.choose.n.cell.type){
  #   sce <- BayesSpace::qTune(sce = sce, qs = seq(n.min, n.max), platform = platform, nrep = nrep, gamma = gamma)
  #   logliks <- attr(sce, "q.logliks")
  #   opt_k <- logliks$q[which(logliks$loglik == max(logliks$loglik))]
  #   sce <- BayesSpace::spatialCluster(sce = sce, q = opt_k, platform = platform, nrep = nrep, gamma = gamma)
  #   #BayesSpace::clusterPlot(sce)
  #   res.cluster <- sce$spatial.cluster
  #   names(res.cluster) <- rownames(sce@colData)
  # }else{
  sce <- BayesSpace::spatialCluster(sce = sce, q = n.setting, d = n.PCs, platform = platform, nrep = nrep, gamma = gamma)
  # BayesSpace::clusterPlot(sce)
  res.cluster <- sce$spatial.cluster
  names(res.cluster) <- rownames(sce@colData)
  # }
  return(res.cluster)
}

#' platform: 'Visium', "ST", "Other_SRT", "scRNAseq"
apply_DR.SC <- function(Seurat.data, n_SVGs = 600, latent_q = 15, num.ct = 7, platform = "Visium"){
  ### extract dataset
  exp_data <-  Seurat.data@assays$SCT@counts
  df <- Seurat.data@meta.data[, c("x", "y", "pixel_x", "pixel_y")]
  temp_colnames <- colnames(df)
  temp_colnames <- sub("x", "row", temp_colnames)
  temp_colnames <- sub("y", "col", temp_colnames)
  colnames(df) <- temp_colnames
  options(Seurat.object.assay.version = "v3")
  Seurat.DR.SC <- Seurat::CreateSeuratObject(counts = exp_data, meta.data = df)
  Seurat.DR.SC <- Seurat::NormalizeData(Seurat.DR.SC, verbose = FALSE)
  Seurat.DR.SC <- DR.SC::FindSVGs(Seurat.DR.SC, nfeatures = n_SVGs, verbose = FALSE)
  Seurat.DR.SC <- DR.SC::DR.SC(seu = Seurat.DR.SC, q = latent_q, K = num.ct, platform = platform)
  # spatialPlotClusters(Seurat.DR.SC)
  res.cluster <- Seurat.DR.SC$spatial.drsc.cluster
  return(res.cluster)
}


apply_SpaGCN <- function(Seurat.data, SpaGCN_beta = 55, SpaGCN_alpha = 1, SpaGCN_p = 0.5,
                         SpaGCN_l = 0.5, num.cluster = 7, res.setting = NULL, SpaGCN_lr = 0.05,
                         SpaGCN_epoches = 200, platform = "Visium"){
  temp <- Seurat.data@images
  counts_path <- temp$counts_path
  meta_path <- temp$loc_path
  img_path <- temp$img_path

  ### for install python scripts in R package
  # EnSDD_pyexec(pyfile = system.file("python", "apply_spaGCN.py",
  #                                   package = "EnSDD"))
  reticulate::source_python(system.file("python", "apply_spaGCN.py",
                                    package = "EnSDD"))
  # reticulate::source_python('/home/vision/Downloads/LHS/EnDecon/EnSDD/inst/python/apply_spaGCN.py')


  run_spaGCN(counts_path, meta_path, img_path, platform = platform, res_setting_by_user = res.setting,
             b = as.integer(SpaGCN_beta), a = SpaGCN_alpha, p = SpaGCN_p,
             n_clusters = as.integer(num.cluster), lr = SpaGCN_lr,
             max_epochs = as.integer(SpaGCN_epoches))

  res.path <- paste0(dirname(img_path), "/spaGCN_res", "/spaGCN_label.txt")
  res_temp <- read.table(res.path, sep = "\t")
  res <- res_temp$V2
  names(res) <- res_temp$V1
  if(sum(res == 0) != 0){
    res <- res + 1
  }
  res <- res[colnames(Seurat.data)]
  unlink(paste0(dirname(img_path), "/spaGCN_res"), recursive=TRUE)
  return(res)
}


apply_stLearn <- function(Seurat.data, n.PCs = 30, n.setting = 7,
                          normalize_type = "physical_distance"){
  temp <- Seurat.data@images
  counts_path <- temp$counts_path
  meta_path <- temp$loc_path
  img_path <- temp$img_path

  ### reticulate source run_stlearn
  # reticulate::source_python('/home/lihs/bio_software/code/apply_stlearn.py')
  # EnSDD_pyexec(pyfile = system.file("python", "apply_stlearn.py",
  #                                   package = "EnSDD"))
  reticulate::source_python(system.file("python", "apply_stlearn.py",
                                        package = "EnSDD"))
  # reticulate::source_python('/home/vision/Downloads/LHS/EnDecon/EnSDD/inst/python/apply_stlearn.py')

  run_stlearn(counts_path, meta_path, img_path, n_comps = as.integer(n.PCs),
              n_setting = as.integer(n.setting),  normalize_type = normalize_type)

  res.path <- paste0(dirname(img_path), "/stLearn_res", "/stLearn_kmeans.txt")
  res_temp <- read.table(res.path, sep = "\t")
  res <- res_temp$V2
  names(res) <- res_temp$V1
  res <- res[colnames(Seurat.data)]
  # }

  unlink(paste0(dirname(img_path), "/stLearn_res"), recursive=TRUE)
  return(res)
}

apply_GraphST <- function(Seurat.data, n_setting = 7,lambda1 = 10, lambda2 = 1,
                          tool = "mclust", radius = 50, n_PCs = 20, refinement = TRUE){
  temp <- Seurat.data@images
  counts_path <- temp$counts_path
  meta_path <- temp$loc_path
  img_path <- temp$img_path

  # reticulate::source_python('/home/lihs/bio_software/code/apply_GraphST.py')
  # EnSDD_pyexec(pyfile = system.file("python", "apply_GraphST.py",
  #                                   package = "EnSDD"))
  reticulate::source_python(system.file("python", "apply_GraphST.py",
                                        package = "EnSDD"))
  # reticulate::source_python('/home/vision/Downloads/LHS/EnDecon/EnSDD/inst/python/apply_GraphST.py')

  if(tool == "mclust"){
    run_GraphST(counts_path, meta_path, img_path, n_setting = as.integer(n_setting),
                lambda1 = lambda1, lambda2 = lambda2, tool = "mclust", radius = as.integer(radius),
                n_PCs = as.integer(n_PCs), refinement = refinement)
    LABEL_PATH <- paste0(dirname(img_path), "/GraphST_res", "/GraphST_mclust.txt")
    pc_temp <- read.table(LABEL_PATH, header = TRUE, sep = "\t")
    sample.names <- pc_temp$X
    pc_temp$X <- NULL
    pc_mat <- as.matrix(pc_temp)
    rownames(pc_mat) <- sample.names
    set.seed(2023)
    library(mclust)
    res <- mclust::Mclust(pc_mat, n_setting, 'EEE')
    res.GraphST <- res$classification
  }else{
    run_GraphST(counts_path, meta_path, img_path, n_setting = as.integer(n_setting),
                lambda1 = lambda1, lambda2 = lambda2, tool = "mclust", radius = as.integer(radius),
                n_PCs = as.integer(n_PCs), refinement = refinement)
    LABEL_PATH <- paste0(dirname(img_path), "/GraphST_res", "/Graph_", tool, ".txt")
    res_temp <- read.table(res.path, sep = "\t")
    res.GraphST <- res_temp$V2
    names(res.GraphST) <- res_temp$V1
    res.GraphST <- res.GraphST[colnames(Seurat.data)]

  }
  unlink(paste0(dirname(img_path), "/GraphST_res"), recursive=TRUE)
  return(res.GraphST)
}

apply_STAGATE <- function(Seurat.data, n_setting = 7, alpha = 0){
  temp <- Seurat.data@images
  counts_path <- temp$counts_path
  meta_path <- temp$loc_path
  img_path <- temp$img_path

  # reticulate::source_python('/home/lihs/bio_software/code/apply_STAGATE.py')
  # EnSDD_pyexec(pyfile = system.file("python", "apply_STAGATE.py",
  #                                   package = "EnSDD"))
  reticulate::source_python(system.file("python", "apply_STAGATE.py",
                                        package = "EnSDD"))
  # reticulate::source_python('/home/vision/Downloads/LHS/EnDecon/EnSDD/inst/python/apply_STAGATE.py')

  run_STAGATE(counts_path, meta_path, img_path, alpha = alpha)
  LABEL_PATH <- paste0(dirname(img_path), "/STAGATE_res", "/STAGATE_mclust.txt")
  pc_temp <- read.table(LABEL_PATH, header = TRUE, sep = "\t")
  sample.names <- pc_temp$X
  pc_temp$X <- NULL
  pc_mat <- as.matrix(pc_temp)
  rownames(pc_mat) <- sample.names
  set.seed(2023)
  library(mclust)
  res <- mclust::Mclust(pc_mat, n_setting, 'EEE')
  res.STAGATE <- res$classification
  unlink(paste0(dirname(img_path), "/STAGATE_res"), recursive=TRUE)
  return(res.STAGATE)
}

# apply_spaVAE <- function(){
#
# }
#
# apply_SiGra <- function(){
#
#
# }

#' Running each base clustering method individually to obtain the  base clustering results on spatially resolved transcriptomics data.
#'
#' This function is implemented to perform individual clustering methods. The current implementation of
#' EnSDD integrates six state-of-the-art methods:  BayesSpace, DR-SC, GraphST, SpaGCN, stLearn and STAGATE.
#' These packages will be automatically installed along with EnSDD. \cr
#'
#' @import reticulate
#' @importFrom Seurat CreateSeuratObject SCTransform FindVariableFeatures ScaleData RunPCA NormalizeData
#' @importFrom abind abind
#' @importFrom Seurat FindNeighbors FindClusters
#' @importFrom SingleCellExperiment SingleCellExperiment
#' @importFrom scran modelGeneVar getTopHVGs
#' @importFrom scater runPCA
#' @importFrom BayesSpace spatialPreprocess qTune spatialCluster
#' @importFrom DR.SC FindSVGs DR.SC
#' @importFrom mclust Mclust
#' @importFrom reticulate source_python use_python py_config
#'
#'
#' @param counts_path the .txt path for SRT gene expression data. The file need include gene names and spots names.
#' @param loc_path the .txt path for meta data of the spots.
#' @param img_path the .tif path of H\& image of the corresponding sequencing tissue.
#' @param python_env the path of python environment. We recommend user construct python environment by the .yml provided by ours.
#' @param n_HVG the number of Highly Variation Genes (HVGs) selected by Seurat. Default setting is 2000.
#' @param n_PCA the number of low-dimensional of SRT gene expression. Default setting is 20.
#' @param number.setting.by.user the number of clustering setting by users. Default setting is 7.
#' @param saving_results whether to save the results of individual clustering results.
#' @param BayesSpace  a logical variable whether to apply BayesSpace.
#' @param DR.SC  a logical variable whether to apply DR.SC.
#' @param SpaGCN  a logical variable whether to apply SpaGCN.
#' @param stLearn a logical variable whether to apply stLearn.
#' @param GraphST  a logical variable whether to apply GraphST.
#' @param STAGATE  a logical variable whether to apply STAGATE.
#' @param BayesSpace_nrep integer variable indicates the number of MCMC iteration.
#' @param DR.SC_n_SVGs integer variable indicates the number of spatially variable genes to be chosen. Default setting is 600.
#' @param DR.SC_latent_q integer variable specifies the number of latent features to be extracted. Default setting is 15.
#' @param SpaGCN_beta integer variable determines the area of each spot when extracting color intensity. Default value is 1.
#' @param SpaGCN_alpha floor variable determines the weight given to histology when calculating Euclidean distance between every two spots. Default setting is 1.
#' @param SpaGCN_p floor variable determines percentage of total expression contributed by neighborhoods.
#' @param SpaGCN_l floor variable control parameter SpaGCN_p.
#' @param SpaGCN_res.setting floor variable determines the resolution for Louvain's algorithm in SpaGCN.
#' @param SpaGCN_lr floor variable determines the learning rate in SpaGCN.
#' @param SpaGCN_epoches integer variable determines the maximum epoches for SpaGCN
#' @param stLearn_n.PCs integer variable determines the number of principal components to compute.
#' @param stLearn_normalize_type item determines the using spatial location (S), tissue morphological feature (M) and gene expression (E) information to normalize data in stLearn. Avariable iterms contains "weights_matrix_all","weights_matrix_pd_gd","weights_matrix_pd_md","weights_matrix_gd_md","gene_expression_correlation","physical_distance","morphological_distance". Default setting is "physical_distance".
#' @param GraphST_lambda1 float variable determines optional weight factor to control the influence of reconstruction loss in mapping matrix learning. Default setting is 10.
#' @param GraphST_lambda2 float variable determines optional weight factor to control the influence of contrastive loss in mapping matrix learning. Default setting is 1.
#' @param GraphST_tool item determines the clustering methods introduced for clustering spots. Avariable tools contain 'leiden', 'louvain', 'mclust'. Default setting is 'mclust'.
#' @param GraphST_radius integer variable determines the number of neighbors considered during refinement. The default is 50.
#' @param GraphST_refinement logical variable refines the predicted labels or not. The default is False.
#' @param GraphST_n_PCs integer variables determines the number PCs in the PCA. Default setting is 20.
#' @param STAGATE_alpha floor variable determines the weight of cell type-aware spatial neighbor network. Default setting is 0.
#'
#' @return a list contains all the results, labels and binary spot-spot similarity matrix, inferred by individual clustering methods and the times of running individual methods. The elements of list is a spots labels vactor, a matrix,  spots * spots and a time vector.
#'
#' @export

run_individual_cluster <- function(counts_path, loc_path, img_path, python_env, n_HVG = 2000, n_PCA = 20,
                                   number.setting.by.user = 7,  saving_results = FALSE,
                                   BayesSpace = TRUE, DR.SC = TRUE, SpaGCN = TRUE, stLearn = TRUE,
                                   GraphST = TRUE, STAGATE = TRUE,
                                   BayesSpace_nrep = 50000,
                                   DR.SC_n_SVGs = 600, DR.SC_latent_q = 15,
                                   SpaGCN_beta = 55, SpaGCN_alpha = 1, SpaGCN_p = 0.5, SpaGCN_l = 0.5,
                                   SpaGCN_res.setting = NULL, SpaGCN_lr = 0.05, SpaGCN_epoches = 200,
                                   stLearn_n.PCs = 30, stLearn_normalize_type = "physical_distance",
                                   GraphST_lambda1 = 10, GraphST_lambda2 = 1, GraphST_tool = "mclust",
                                   GraphST_radius = 50, GraphST_refinement = TRUE, GraphST_n_PCs = 20,
                                   STAGATE_alpha = 0){
  if(is.null(python_env)){
    cat("if you want to use the environment of your own python, please set it manual.
         Otherwise, the code will create R_reticulate environment and install the essential
         python packages for running SpaGCN, stLearn, GraphST, STAGATE.
         We strongly recommend the user to construct anaconda env and install python
         packages by running the .yml file we provided.")
  }else{
    reticulate::use_python(python_env, require = T)
    reticulate::py_config()
  }
  message("Some methods, including STAGATE, DR.SC, BayesSpace
          needs the number of clusters setting by users, Other methods select the number
          of cluster by various strategies.\n")
  Methods <- c("BayesSpace", "DR.SC", "SpaGCN", "stLearn", "GraphST", "STAGATE")
  Methods.idx <- c(BayesSpace, DR.SC, SpaGCN, stLearn, GraphST, STAGATE)
  Methods.used <- Methods[Methods.idx]
  K <- length(Methods.used)

  Seurat.data <- data_process(counts_path, loc_path, img_path, n_HVG, n_PCA)
  ### global setting
  sequencing_platform = "Visium"

  res.clustering <- list()
  time.methods <- c()
  k <- 1
  # apply Giotto
  # if(Giotto){
  #   cat("Run Giotto_LD, Giotto_KM, Giotto_H, Giotto_HM...\n")
  #   pred.Giotto <- Sys.time()
  #   temp.Giotto <- apply_Giotto(Seurat.data, LD_k = Giotto_LD_k, LD_resolution = Giotto_LD_resolution,
  #                               num_cluster = number.setting.by.user, HM_walk_step = Giotto_HM_walk_step)
  #   end.Giotto <- Sys.time()
  #   time.Giotto <- difftime(end.Giotto, pred.Giotto, units = "mins")
  #   cat("Run time for Giotto: ", time.Giotto, "min","\n")
  #   time.methods <- c(time.methods, time.Giotto)

  #   res.clustering$Giotto_LD <- temp.Giotto$Giotto_LD
  #   res.clustering$Giotto_KM <- temp.Giotto$Giotto_KM
  #   res.clustering$Giotto_H <- temp.Giotto$Giotto_H
  #   res.clustering$Giotto_HM <- temp.Giotto$Giotto_HM
  #   k <- k + 1
  # }

  # apply Seurat
  # if(Seurat){
  #   cat("Run Seurat...\n")
  #   pred.Seurat <- Sys.time()
  #   temp.Seurat <- apply_Seurat(Seurat.data, Seurat_k_neighbor, Seurat.resolution, number.setting.by.user)
  #   end.Seurat <- Sys.time()
  #   time.Seurat <- difftime(end.Seurat, pred.Seurat, units = "mins")
  #   cat("Run time for Seurat: ", time.Seurat, "min","\n")
  #   time.methods <- c(time.methods, time.Seurat)
  #   res.clustering$Seurat <- temp.Seurat
  #   k <- k + 1
  # }
  # apply BayesSpace
  if(BayesSpace){
    cat("Run BayesSpace...\n")
    pred.BayesSpace <- Sys.time()
    temp.BayesSpace <- apply_BayesSpace(Seurat.data, platform = sequencing_platform, n.PCs = n_PCA,
                                        n.HVGs = n_HVG, n.setting = number.setting.by.user,
                                        nrep = BayesSpace_nrep)
    end.BayesSpace <- Sys.time()
    time.BayesSpace <- difftime(end.BayesSpace, pred.BayesSpace, units = "mins")
    cat("Run time for BayesSpace: ", time.BayesSpace, "min","\n")
    time.methods <- c(time.methods, time.BayesSpace)
    res.clustering$BayesSpace <- temp.BayesSpace
    k <- k + 1
  }
  # apply DR.SC
  if(DR.SC){
    cat("Run DR.SC...\n")
    pred.DR.SC <- Sys.time()
    temp.DR.SC <- apply_DR.SC(Seurat.data, n_SVGs = DR.SC_n_SVGs, latent_q = DR.SC_latent_q,
                              num.ct = number.setting.by.user, platform = sequencing_platform)
    end.DR.SC <- Sys.time()
    time.DR.SC <- difftime(end.DR.SC, pred.DR.SC, units = "mins")
    cat("Run time for DR.SC: ", time.DR.SC, "min","\n")
    time.methods <- c(time.methods, time.DR.SC)
    res.clustering$DR.SC <- temp.DR.SC
    k <- k + 1
  }
  # apply SpaGCN
  if(SpaGCN){
    cat("Run SpaGCN...\n")
    pred.SpaGCN <- Sys.time()
    temp.SpaGCN <- apply_SpaGCN(Seurat.data, SpaGCN_beta, SpaGCN_alpha, SpaGCN_p,
                                SpaGCN_l, number.setting.by.user, SpaGCN_res.setting, SpaGCN_lr,
                                SpaGCN_epoches, sequencing_platform)
    end.SpaGCN <- Sys.time()
    time.SpaGCN <- difftime(end.SpaGCN, pred.SpaGCN, units = "mins")
    cat("Run time for SpaGCN: ", time.SpaGCN, "min","\n")
    time.methods <- c(time.methods, time.SpaGCN)
    res.clustering$SpaGCN <- temp.SpaGCN
    k <- k + 1
  }
  # apply stLearn
  if(stLearn){
    cat("Run stLearn...\n")
    pred.stLearn <- Sys.time()
    temp.stLearn <- apply_stLearn(Seurat.data, stLearn_n.PCs, number.setting.by.user, stLearn_normalize_type)
    end.stLearn <- Sys.time()
    time.stLearn <- difftime(end.stLearn, pred.stLearn, units = "mins")
    cat("Run time for stLearn: ", time.stLearn, "min","\n")
    time.methods <- c(time.methods, time.stLearn)
    res.clustering$stLearn <- temp.stLearn
    k <- k + 1
  }
  # apply GraphST
  if(GraphST){
    cat("Run GraphST...\n")
    pred.GraphST <- Sys.time()
    temp.GraphST <- apply_GraphST(Seurat.data, number.setting.by.user, GraphST_lambda1, GraphST_lambda2,
                                  GraphST_tool, GraphST_radius, GraphST_n_PCs, GraphST_refinement)
    end.GraphST <- Sys.time()
    time.GraphST <- difftime(end.GraphST, pred.GraphST, units = "mins")
    cat("Run time for GraphST: ", time.GraphST, "min","\n")
    time.methods <- c(time.methods, time.GraphST)
    res.clustering$GraphST <- temp.GraphST
    k <- k + 1
  }
  # apply STAGATE
  if(STAGATE){
    cat("Run STAGATE...\n")
    pred.STAGATE <- Sys.time()
    temp.STAGATE <- apply_STAGATE(Seurat.data, number.setting.by.user, STAGATE_alpha)
    end.STAGATE <- Sys.time()
    time.STAGATE <- difftime(end.STAGATE, pred.STAGATE, units = "mins")
    cat("Run time for STAGATE: ", time.STAGATE, "min","\n")
    time.methods <- c(time.methods, time.STAGATE)
    res.clustering$STAGATE <- temp.STAGATE
    k <- k + 1
  }

  names(time.methods) <- Methods.used
  ### check NA in the methods
  index.list <- c()
  for (i in 1:length(res.clustering)) {
    na.index <- sum(is.na(res.clustering[[i]]))
    if(na.index != 0){
      index.list <- c(index.list, i)
    }
  }

  if(length(index.list) != 0){
    methods.names.tmp <- names(res.clustering)[index.list]
    for (i in methods.names.tmp) {
      res.clustering[[i]] <- NULL
    }
  }
  ### check NA in the methods
  res.clustering.matrix <- lapply(res.clustering, clust2Mat)
  res <- list(clustering.vec = res.clustering , clustering.mat = res.clustering.matrix, time = time.methods)
  if(saving_results)
  {
    cat("Saving individual clustering analysis results...\n")

    save(res, file = "res.clustering.RData")
  }
  return(res)
}
