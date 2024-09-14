


DotPlot <- function(x, #SeuratObj; dgCMatrix; Matrix
                    assay, #input assay/ output of searchMarker
                    genes,
                    dot.scale){
  UseMethod("DotPlot")
}

#' @import Seurat
#' @import magrittr
#' @import dplyr
#' @export
DotPlot.Seurat <- function(x, #SeuratObj; dgCMatrix; Matrix
                           assay, #input assay/ output of searchMarker
                           genes,
                           dot.scale){
  #core calculation
  p <- DotPlot(
    x,
    assay = "RNA",
    features = genes,
    #cols = c("white", "red"), #replaced with scale color function
    col.min = -2.5,
    col.max = 2.5,
    dot.min = 0,
    dot.scale = dot.scale,
    idents = NULL,
    group.by = NULL,
    split.by = NULL,
    cluster.idents = FALSE,
    scale = TRUE,
    scale.by = "radius",
    scale.min = NA,
    scale.max = NA
  ) + coord_flip() + scale_colour_gradient2(low = colors.gradient[1],mid = colors.gradient[2], high = colors.gradient[3]) +
    ggthemes::theme_few() + theme(axis.text = element_text(size = 8),axis.ticks = element_blank(),axis.text.x = element_text(angle = 45, hjust = 1))
}
