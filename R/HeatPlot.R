#' HeatPlot
#'
#' @param x a seurat object
#' @param assay the assay in the Seurat Object to use
#' @param genes the genes for the Seurat Object to use
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
                     genes){

  mat <- Seurat::AverageExpression(x,assays = assay)[[1]][genes,]
  mat <- log10(mat + 1)

  p <- pheatmap::pheatmap(mat,color = viridis::viridis(10),border_color = NA)

  return(p)
}
