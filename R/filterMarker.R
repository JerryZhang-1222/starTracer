#' filterMarker
#'
#' @param x Seurat Object
#' @param ident.use ident to use during calculation, should be the same as the one user used during FindAllMarkers()
#' @param mat the out out data.frame of FindAllMarkers()
#' @param num the number of marker genes the user prefers
#' @param thresh.min the minimum value that consider a given cell expresses a certain gene
#'
#' @import tibble
#' @import tidyr
#' @import dplyr
#' @import magrittr
#' @import rlang
#'
#' @return a modified data.frame of the FindAllMarkers function
#'
#' @export
#'
#' @examples \dontrun{filterMarker(x = pbmc_small, ident.use = "RNA_snn_res.1", mat = diff.wilcox)}
filterMarker <- function(x,ident.use,mat,num = "all",thresh.min = 0){
  #get num
  if(num == "all"){
    message("num is set to all, now finding the optimal number...")
    num <- max(as.matrix(table(mat$cluster))[,1])
    message(paste0("using ",num," as the maximumn number to find marker genes in each cluster"))
  }


  #calculate pct.pos
  if(ncol(x) < 100000){
    mat <- cal_pctpos(x = x, ident.use = ident.use, mat = mat, thresh.min = 0)
  } else {
    message("detecting cells greater than 100,000, using pct.pos calculating method for big data...")
    mat <- cal_pctpos.big(x = x, ident.use = ident.use, mat = mat, thresh.min = 0)
  }

  mat <- mat[order(-mat[, "pct.pos"]),]

  mat <- mat[!duplicated(mat$gene),] %>% arrange(cluster,desc(pct.pos),desc(avg_log2FC))

  mat <- mat %>%
    group_by(cluster) %>%
    slice_head(n = num) %>%
    ungroup()

  #message:
  info <- mat$cluster %>% table() %>% as.data.frame()
  info$echo <- dplyr::case_when(info$Freq == num ~ "Yes",
                                info$Freq != num ~ "NO")
  info.yes <- paste(subset(info,info$echo == "Yes")[,1],collapse = ",")
  info.no <- paste(subset(info,info$echo == "NO")[,1],collapse = ",")
  info.no.num <- paste(subset(info,info$echo == "NO")[,2],collapse = ",")
  if(info.yes == ""){
    warning(paste0("We found no cluster with ", num, " marker genes in each, please make sure you have the correct input data or try to assign a smaller num"))
  } else {
    message(paste0("We found ",num," marker genes for cluster: ",info.yes))
  }
  if(info.no != ""){
    message(paste0("starTracer found ", num," marker genes in cluster: ",info.yes,", for clusters ", info.no, " starTracer found ", info.no.num, " marker genes."))
  }

  genes.markers <- mat$gene
  p <- HeatPlot(x = x, assay = "RNA", genes = genes.markers)

  filtermarkers.out <- list(marker.frame = mat, heatmap = p, genes.markers = genes.markers)

  return(filtermarkers.out)
}
