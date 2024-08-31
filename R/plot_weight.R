library(ggpubr)
library(ggrepel)

#' plot weight and ARI dotplot
#'
#' @importFrom ggplot2 ggplot ggsave
#' @importFrom ggrepel geom_text_repel
#' @importFrom ggpubr stat_cor
#' @param data_bar the dataframe for plotting of dotplot.
#' @param sample_id the title for the dotplot.
#' @param file_name the file name of the saving plot.
#' @param save_file whether save the dotplot.
#' @param label.x the x-axis location of the legend.
#' @param label.y.inc the y-axis location of the legend.


plot_weight <- function(data_bar, sample_id, file_name, save_file = FALSE,
                        label.x = 0.12, label.y.inc = 0.03){
  title_sample <-  sample_id
  p2 = ggplot2::ggplot(data_bar, aes(y = ARI, x = weight_est, color = factor(Method))) +
    geom_point(size = 3) +
    theme_classic() +
    labs(x="Weights assigned to individual methods",
         y="ARI of individual methods",title = title_sample) +
    theme(axis.text.x = element_text(size = 12, hjust = 1),
          axis.text.y.left = element_text(size = 12),
          axis.title.y = element_text(size = 14),
          axis.title.x = element_text(size = 14),
          panel.border = element_blank(),
          legend.position = "none",
          axis.text.y = element_text(size = 14),
          plot.title = element_text(hjust = 0.5, size = 14)) +
    ggrepel::geom_text_repel(aes(label = Method), size = 5) +
    # stat_cor(method='pearson', cor.coef.name = "R", size = 4,
    #          label.y = min(data_bar$ARI) + 0.15, label.x = 0.12,color='red') +
    # stat_cor(method='spearman', cor.coef.name = 'rho', size = 4,
    #          label.y = min(data_bar$ARI) + 0.1, label.x = .12,color='red')
    ggpubr::stat_cor(method='pearson', cor.coef.name = "R", size = 5,
                     label.y = min(data_bar$ARI) + label.y.inc + 0.04, label.x = label.x,color='red') +
    stat_cor(method='spearman', cor.coef.name = 'rho', size = 5,
             label.y = min(data_bar$ARI) + label.y.inc, label.x = label.x,color='red')
  if(save_file){
    ggplot2::ggsave(p2, filename = file_name, width = 6, height = 5)
  }else{
    return(p2)
  }
  
}

