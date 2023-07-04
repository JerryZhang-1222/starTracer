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
  meta <- FetchData(x, c(ident.use)) %>%
    rownames_to_column("cell.id")
  data <- GetAssayData(x, slot = "counts", assay = "RNA")
  ident_sym <- sym(ident.use)
  mat[,"pct.pos"] <- NA

  for(i in 1:nrow(mat)) {
    ident <- mat[i,"cluster"] %>% as.character()
    gene <- mat[i,"gene"]
    line <- data[gene,, drop = F]
    expr.all <- rowSums(line > thresh.min) #5846
    expr.id <- rowSums(line[,which(colnames(line) %in% meta[meta[,ident.use] == ident,"cell.id"]),drop = F] > thresh.min)
    mat[i,"pct.pos"] <- expr.id/expr.all
  }
  return(mat)
}

