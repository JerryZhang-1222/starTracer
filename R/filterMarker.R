#' filterMarker
#'
#' @param x Seurat Object
#' @param ident.use ident to use during calculation, should be the same as the one user used during FindAllMarkers()
#' @param mat the out out data.frame of FindAllMarkers()
#'
#' @return a modified data.frame of the FindAllMarkers function
#' @export
#'
#' @examples \dontrun{filterMarker(x = pbmc_small, ident.use = "RNA_snn_res.1", mat = diff.wilcox)}
filterMarker <- function(x,ident.use,mat){

  ref <- x@meta.data[,ident.use] %>% table() %>% as.data.frame()
  colnames(ref) <- c("ident_use","count")
  rownames(ref) <- ref[,"ident_use"]
  ref[,"ident_use"] <- NULL

  order <- c()
  for(i in 1:nrow(mat)){
    line <- mat[i,]
    ident <- line["cluster"] %>% as.character()
    ident.num <- ref[ident,"count"] %>% as.numeric()
    rest.num <- as.numeric(sum(ref$count) - ident.num)
    pct1 <- line["pct.1"] %>% as.numeric()
    pct2 <- line["pct.2"] %>% as.numeric()
    rate.i <- (ident.num*pct1)/(ident.num*pct1 + rest.num*pct2)
    order <- append(order,rate.i)
  }

  mat$trace_order <- order
  mat <- mat[order(mat[,"cluster"],-mat[,"trace_order"]),]

  return(mat)
}
