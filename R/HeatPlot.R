#' HeatPlot
#'
#' @param x a seurat object
#' @param assay the assay in the Seurat Object to use
#' @param genes the genes for the Seurat Object to use
#' @param scale.method the method use to scale the original matrix
#' @param lim.scale when using "scale" as the method of scale.method, scaled value will be limited from -lim.scale to lim.scale
#' @param border_color the color of border in the heatmap
#'
#' @return a plot of heatmap
#'
#' @import Seurat
#' @import pheatmap
#' @import viridis
#'
#' @examples \dontrun{HeatPlot(x = pbmc_small,assay = "RNA",genes = genes)}
HeatPlot <- function(x,
                     assay,
                     genes,
                     scale.method = "scale",
                     lim.scale = 2,
                     border_color = "black"){
  mat <- Seurat::AverageExpression(x,assays = assay)[[1]][genes,]
  if(scale.method=="log"){
    mat <- log10(mat + 1)
  } else if(scale.method=="scale"){
    mat <- t(mat) %>% scale() %>% t()
    mat[mat > lim.scale] <- lim.scale
    mat[mat < -lim.scale] <- -lim.scale
  }
  p <- pheatmap::pheatmap(mat,color = viridis::viridis(10),border_color = border_color)

  return(p)
}
