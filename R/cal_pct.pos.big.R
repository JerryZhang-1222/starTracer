#' cal_pctpos.big
#'
#' @param x the seurat object
#' @param ident.use the ident to use
#' @param mat marker gene table of seurat
#' @param thresh.min
#'
#' @import tibble
#' @import Seurat
#' @import rlang
#' @import dplyr
#' @import Matrix
#'
#' @return a modified marker gene matrix with a new colum named "pct.pos"
#'
#' @examples \dontrun{cal_pctpos.big(x = x,ident.use = ident.use,mat = mat,thresh.min = 0)}
cal_pctpos.big <- function(x,ident.use,mat,thresh.min = 0){
  message("Fetching meta.data")
  meta <- FetchData(x, c(ident.use)) %>%
    rownames_to_column("cell.id")

  data <- GetAssayData(x, slot = "counts", assay = "RNA")
  message("converting row to column...")
  data.t <- Matrix::t(data)

  t1 <- Sys.time()
  data.t[,1,drop = F]
  t2 <- Sys.time()
  t <- t2 - t1
  est.time <- round(as.numeric(t)*nrow(mat)/60,digits = 3)

  message(paste0("ESTIMATING TIME LEFT: ",est.time," min"))
  if(est.time > 5){message("why not take a rest and have a cup of tea?")}

  mat[,"pct.pos"] <- NA

  message("now calculating pct.pos...")
  for(i in 1:nrow(mat)) {
    if(i %% 100 == 0){message(paste0("calculating pct.pos for the ",i,"th gene..."))}
    ident <- mat[i,"cluster"] %>% as.character()
    gene <- mat[i,"gene"]
    line <- data.t[,gene, drop = F]
    expr.all <- colSums(line > thresh.min) #5846
    expr.id <- colSums(line[which(rownames(line) %in% meta[meta[,ident.use] == ident,"cell.id"]),,drop = F] > thresh.min)
    mat[i,"pct.pos"] <- expr.id/expr.all
  }
  return(mat)
}

