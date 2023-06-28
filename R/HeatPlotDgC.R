#' HeatPlotDgC
#'
#' @param x the average expression matrix
#' @param genes genes we will use
#'
#' @return a heatmap
#'
#' @examples \dontrun{HeatPlotDgC(x = pbmc_small,genes = genes)}
HeatPlotDgC <- function(x, genes){
  mat <- x[genes,]
  mat <- t(mat) %>% scale() %>% t()
  pheatmap::pheatmap(mat)
}
