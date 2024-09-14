#' filtGene_nagative
#'
#' this function is used to filter genes in each cluster acoording to the given thresh.2.neg when calculating negative markers
#' @param expr.norm.neg the max.nomalized expression matrix for negative marker
#' @param thresh.2.neg this is the threshold range from 0 to 1. thresh.2.neg will be automatically mapped to the different quantiles of the range of the mean expression level in each different cluster
#' @param clusters a vector of all clusters in the single cell data
#'
#' @return a matrix fater filtering
#'
#' @import magrittr
#' @import dplyr
#'
#' @examples \dontrun{filtGene(expr.norm.neg = x,thresh.2.neg = 0.5)}
filtGene.Neg <- function(expr.norm.neg,thresh.2.neg,clusters){

  #load the max.normalized data:
  mat <- data.frame()

  #for each cluster, keep lines with mean.orig bigger than thresh.2.neg
  for(i in clusters){
    mat.i <- subset(expr.norm.neg,expr.norm.neg[,"max.X"] == i)
    mat.i <- arrange(mat.i,mat.i[,"mean.orig"])

    mat.i <- mat.i[1:round(nrow(mat.i)*(1-thresh.2.neg)),]

    mat <- rbind(mat,mat.i)
  }

  return(mat)
}
