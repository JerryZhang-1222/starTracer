#' filterMarker
#'
#' @param x Seurat Object
#' @param ident.use ident to use during calculation, should be the same as the one user used during FindAllMarkers()
#' @param mat the out out data.frame of FindAllMarkers()
#' @param num the number of marker genes the user prefers
#'
#' @return a modified data.frame of the FindAllMarkers function
#' @export
#'
#' @examples \dontrun{filterMarker(x = pbmc_small, ident.use = "RNA_snn_res.1", mat = diff.wilcox)}
filterMarker <- function(x,ident.use,mat,num = "all"){
  #get num
  if(num == "all"){
    message("num is set to all, now finding the optimal number...")
    num <- max(as.matrix(table(mat$cluster))[,1])
    message(paste0("using ",num," as the maximumn number to find marker genes in each cluster"))
  }
  #
  ref <- x@meta.data[, ident.use] %>% table() %>% as.data.frame()
  colnames(ref) <- c("ident_use", "count")
  rownames(ref) <- ref[, "ident_use"]
  ref[, "ident_use"] <- NULL
  order <- c()
  for (i in 1:nrow(mat)) {
    line <- mat[i, ,drop = T]
    ident <- line$cluster %>% as.character()
    ident.num <- ref[ident, "count"] %>% as.numeric()
    rest.num <- as.numeric(sum(ref$count) - ident.num)
    pct1 <- line["pct.1"] %>% as.numeric()
    pct2 <- line["pct.2"] %>% as.numeric()
    rate.i <- (ident.num * pct1)/(ident.num * pct1 + rest.num *
                                    pct2)
    order <- append(order, rate.i)
  }
  mat$pct.pos <- order
  mat <- mat[order(-mat[, "pct.pos"]),]

  mat <- mat[!duplicated(mat$gene),] %>% arrange(cluster,desc(pct.pos))

  mat <- mat %>%
    group_by(cluster) %>% #按照cluster列进行分组
    top_n(num, wt = pct.pos) %>% #按照FC列从大到小排序，选择每组前2行
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
