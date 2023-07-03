#' cal_pctpos
#'
#' @param x the seurat object
#' @param ident.use the ident to use
#' @param mat marker gene table of seurat
#' @param thresh.min
#'
#' @return a modified marker gene matrix with a new colum named "pct.pos"
#'
#' @examples \dontrun{cal_pctpos(x = x,ident.use = ident.use,mat = mat,thresh.min = 0)}
cal_pctpos <- function(x,ident.use,mat,thresh.min = 0){
  meta <- FetchData(x, c(ident.use)) %>%
    tibble::rownames_to_column("cell.id")
  data <- GetAssayData(x, slot = "counts", assay = "RNA")

  .tmp <- as.data.frame(t(as.matrix(example.SO.counts[which(rownames(data) %in% mat$gene),])))
  .tmp <- tibble::rownames_to_column(.tmp,"cell.id")
  .tmp <- left_join(.tmp,meta, by = "cell.id")
  .tmp <- gather(.tmp, key = "gene", value = "nUMI", c(-cell.id, -cell_type__ontology_label))
  .tmp <- .tmp %>% group_by(cell_type__ontology_label, gene) %>%
    dplyr::summarise(non.expressing.cells = sum(nUMI == 0), expressing.cells = sum(nUMI > 0)) %>%
    ungroup()

  .tmp <- subset(.tmp,paste(.tmp$gene,.tmp$cell_type__ontology_label,sep = "_") %in% paste(mat$gene,mat$cluster,sep = "_"))
  .tmp$expressing.cells.all <- rowSums(data[.tmp$gene, ,drop = FALSE] > thresh.min)

  .tmp <- dplyr::mutate(.tmp,pct.pos = expressing.cells/expressing.cells.all)
  .tmp <- .tmp[match(mat$gene,.tmp$gene),]

  mat$pct.pos <- .tmp$pct.pos
  return(mat)
}

