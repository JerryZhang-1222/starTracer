#' StVlnPlot
#'
#' @param x a seurat object
#' @param assay the assay to use, set as "RNA" by default
#' @param genes genes to draw plot, we recommand to use the output gene list of starTracer::searchMarker()
#' @param ident the ident to use to draw a plot, the names must be a column from the seuratobject@@metadata
#' @param colors colors to use
#'
#' @return a stackedVionlnPlot
#' @import dplyr
#' @import reshape2
#' @import ggplot2
#' @import viridis
#' @import magrittr
#'
#'
#' @export
#'
#' @examples
#' \dontrun{
#' StVlnPlot(x = pbmc_small,
#' assay = "RNA",
#' genes = genes,
#' ident = "RNA_snn_res.1",
#' colors = NULL)
#' }
StVlnPlot <- function(x,assay = "RNA",genes,ident = "seurat_clusters",colors = NULL){
  vln.df <- as.data.frame(x[[assay]]@data[genes,])
  vln.df$gene <- rownames(vln.df)

  vln.df <- reshape2::melt(vln.df,id="gene")
  colnames(vln.df)[c(2,3)] <- c("cells","exp")

  anno <- x@meta.data[,c("cells",ident)]

  vln.df <- dplyr::inner_join(vln.df,anno,by="cells")
  vln.df$gene <- factor(vln.df$gene,levels = genes)

  n <- x@meta.data[,ident] %>% unique() %>% length()
  if(is.null(colors)){colors <- viridis::viridis(n)}
  #plot
  p <- vln.df %>% ggplot(aes_string("exp",ident)) +
    ggplot2::geom_violin(aes_string(fill=ident), scale = "width") +
    ggplot2::facet_grid(.~gene, scales = "free_x") +
    ggplot2::scale_fill_manual(values = colors) +
    ggplot2::scale_x_continuous(expand = c(0,0)) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.title.y.left = ggplot2::element_blank(),
      axis.ticks.y.left = ggplot2::element_blank(),
      axis.text.y.left = ggplot2::element_text(angle = 0,hjust = 1,vjust = 0.5,color = "black",size = 14),
      axis.title.x.bottom = ggplot2::element_blank(),
      axis.ticks.x.bottom = ggplot2::element_blank(),
      axis.text.x.bottom = ggplot2::element_blank(),
      legend.position = "none",
      panel.spacing.x = unit(0, "cm"),
      strip.text.x = ggplot2::element_text(angle= 90,size = 14,hjust = 0),
      strip.background.x = ggplot2::element_blank()
    )

  return(p)
}
