# data = read.table("C:/data/Vandy Onedrive/OneDrive - Vanderbilt/Desktop/Coranl_df.txt")

create_df_Sankey <- function(cluster1, cluster2){
  c1toc2 <- NULL

  levels1 <- unique(cluster1)
  levels2 <- unique(cluster2)

  for(i in 1:length(levels1)){
    res <- NULL
    for(j in 1:length(levels2)) {
      tmp <- length(which(cluster1 == levels1[i] & cluster2 == levels2[j]))
      res <- c(res, tmp)
    }
    c1toc2 <- rbind(c1toc2, res)
  }
  colnames(c1toc2) <- levels2
  rownames(c1toc2) <- levels1

  return(c1toc2)
}

plot_sankey_fun = function(links){
  nodes <- data.frame(
    name=c(as.character(links$source),
           as.character(links$target)) %>% unique()
  )

  links$IDsource <- match(links$source, nodes$name)-1
  links$IDtarget <- match(links$target, nodes$name)-1

  p <- networkD3::sankeyNetwork(Links = links, Nodes = nodes,
                                Source = "IDsource", Target = "IDtarget",
                                Value = "value", NodeID = "name",
                                sinksRight=FALSE, fontSize = 12)
}

#' This function is designed for the visualization of Sankey diagram.
#'
#'
#' @importFrom networkD3 sankeyNetwork
#' @importFrom reshape2 melt
#'
#' @param df the data frame, where row represents the sample, column represents the multiple clusters
#' @param filter filter the number of interaction in Sankey diagram. Default setting is 5.
#'
#'
#' @export
#'

plot_sankey <- function(df, filter = 5){
  data <- df
  df_list <- list()
  # source_sk <- c()
  # target_sk <- c()
  # values_sk <- c()
  for (i in 1:(ncol(df)-1)) {
    tmp1 <- create_df_Sankey(df[,i], df[,(i+1)])

    level_row <- paste0('Domain', length(levels(as.factor(df[,i]))), '_', rownames(tmp1))
    level_col <- paste0('Domain', length(levels(as.factor(df[,i+1]))), '_', colnames(tmp1))

    rownames(tmp1) <- level_row
    colnames(tmp1) <- level_col

    tmp2 <- reshape2::melt(tmp1)
    tmp2 = tmp2[tmp2$value >= filter, ]
    # source_sk <- c(source_sk, tmp2[,1])
    # target_sk <- c(target_sk, tmp2[,2])
    # values_sk <- c(values_sk, tmp2[,3])

    df_list[[i]] <- tmp2
  }

  source_sk <- c()
  target_sk <- c()
  values_sk <- c()
  for (j in 1:length(df_list)) {
    source_sk <- c(source_sk, as.character(df_list[[j]][,1]))
    target_sk <- c(target_sk, as.character(df_list[[j]][,2]))
    values_sk <- c(values_sk, df_list[[j]][,3])
  }

  links <- data.frame(
    source=source_sk,
    target=target_sk,
    value=values_sk
  )
  p = plot_sankey_fun(links)
  p
}

