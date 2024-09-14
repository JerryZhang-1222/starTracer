#' filtGene
#'
#' this function is used to filter genes in each cluster acoording to the given thresh.2
#' @param expr.norm the max.nomalized expression matrix
#' @param thresh.2 this is the threshold range from 0 to 1. Thresh.2 will be automatically mapped to the different quantiles of the range of the mean expression level in each different cluster
#' @param clusters a vector of all clusters in the single cell data
#'
#' @return a matrix fater filtering
#'
#' @import magrittr
#' @import dplyr
#'
#' @examples \dontrun{filtGene(expr.norm = x,thresh.2 = 0.5)}
filtGene <- function(expr.norm,thresh.2,clusters){

  #load the max.normalized data:
  mat <- data.frame()

  #for each cluster, keep lines with mean.orig bigger than thresh.2
  for(i in clusters){
    mat.i <- subset(expr.norm,expr.norm[,"max.X"] == i)
    if(nrow(mat.i) != 0){
      mat.i <- arrange(mat.i,desc(mat.i[,"mean.orig"]))

      if(nrow(mat.i) != 1){mat.i <- mat.i[1:round(nrow(mat.i)*(1-thresh.2)),]}

      mat <- rbind(mat,mat.i)
    }
  }

  return(mat)
}
