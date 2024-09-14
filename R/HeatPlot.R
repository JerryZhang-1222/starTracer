#' HeatPlot
#'
#' @param x a seurat object
#' @param genes the genes for the Seurat Object to use
#' @param scale.method the method use to scale the original matrix
#' @param lim.scale when using "scale" as the method of scale.method, scaled value will be limited from -lim.scale to lim.scale
#' @param border_color the color of border in the heatmap
#' @param colors the colorkey of the heatmap
#' @param clusters the clusters of import data
#'
#' @return a plot of heatmap
#'
#' @import Seurat
#' @import ComplexHeatmap
#' @import viridis
#' @import grid
#'
#' @examples \dontrun{HeatPlot(x = pbmc_small,assay = "RNA",genes = genes)}

HeatPlot <- function(x,
                     genes,
                     genes.top1,
                     scale.method = "scale",
                     lim.scale = 2,
                     border_color = "black",
                     clusters = clusters,
                     colors = viridis(10)){
  mat <- x[genes,]
  if(scale.method=="log"){
    mat <- log10(mat + 1)
  } else if(scale.method=="scale"){
    mat <- t(mat) %>% scale()
    mat[mat > lim.scale] <- lim.scale
    mat[mat < -lim.scale] <- -lim.scale
  }
  rownames(mat) <- clusters
 #plotting heatmap
  marked_indices <- which(colnames(mat) %in% genes.top1)
  p <- Heatmap(mat,
               name = "Expression",
               cluster_rows = F,
               cluster_columns = F,
               show_column_names = FALSE,
               show_row_names = T,
               row_names_gp = gpar(fontsize = 5),
               column_dend_height = unit(4, "cm"),
               heatmap_legend_param = list(title = "z-score"),
               top_annotation = HeatmapAnnotation(
                 mark = anno_mark(
                   at = marked_indices,
                   labels = colnames(mat)[marked_indices]
                 ))
  )
}






























































































