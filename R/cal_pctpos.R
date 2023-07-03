#' cal_pctpos
#'
#' @param x the seurat object
#' @param ident.use the ident to use
#' @param mat marker gene table of seurat
#' @param thresh.min
#'
#' @import tibble
#' @import tidyr
#' @import dplyr
#' @import magrittr
#' @import rlang
#'
#' @return a modified marker gene matrix with a new colum named "pct.pos"
#'
#' @examples \dontrun{cal_pctpos(x = x,ident.use = ident.use,mat = mat,thresh.min = 0)}
cal_pctpos <- function(x,ident.use,mat,thresh.min = 0){
  meta <- FetchData(x, c(ident.use)) %>%
    tibble::rownames_to_column("cell.id")
  data <- GetAssayData(x, slot = "counts", assay = "RNA")

  ident_sym <- sym(ident.use)

  message("getting dgCMatrix from object")
  .tmp <- as.data.frame(t(as.matrix(data[which(rownames(data) %in% mat$gene),])))
  .tmp <- rownames_to_column(.tmp,"cell.id")
  .tmp <- left_join(.tmp,meta, by = "cell.id")

  message("calculating gathering result")
  .tmp <- gather(.tmp, key = "gene", value = "nUMI", c(-cell.id, -!!ident_sym))

  message("calculating expression info")
  .tmp <- .tmp %>% group_by(!!ident_sym, gene) %>%
    dplyr::summarise(non.expressing.cells = sum(nUMI == 0), expressing.cells = sum(nUMI > 0)) %>%
    ungroup()

  .ref <- as.data.frame(rowSums(data[, ,drop = FALSE] > thresh.min)) %>% rownames_to_column("gene")
  colnames(.ref)[2] <- "expressing.cells.all"

  .tmp <- merge(.tmp,.ref,by = "gene")

  .tmp <- mutate(.tmp,pct.pos = expressing.cells/expressing.cells.all)
  .tmp <- .tmp[match(mat$gene,.tmp$gene),]

  mat$pct.pos <- round(.tmp$pct.pos,digits = 3)
  return(mat)
}

