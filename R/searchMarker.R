#' searchMarker
#'
#' this is the main function of searching marker genes
#'
#' @param x a Seurat Object, an average expression matrix or a sparse matrix.
#' @param thresh.1 the threshold for dividing clusters max-normalized gene expression to high or low, set to 0.5 by default
#' @param thresh.2 the specificity level or the lowest expression level of the genes user prefers
#' @param method the arranging method of ordering the genes. Set default and strongly suggestted as "del_MI"
#' @param num the topN gene user prefers
#' @param gene.use the genes we use, set NULL to use all the genes or set as "HVG" to only use HVG in Seurat Object
#' @param meta.data when importing a dgCMatrix and a meta.data, users should pass the meta.data to this param
#' @param ident.use when importing a dgCMatrix and a meta.data, users should specify the ident that should be used in the meta.data
#' @param scale.method the scaling method that will be used when outputting the heatmap
#' @param lim.scale when using "scale" as the method of scale.method, scaled value will be limited from -lim.scale to lim.scale
#' @param border_color the color of border in the heatmap
#' @param colors the colorkey of the heatmap
#'
#' @return a list containing a pseudobulk expression matrix of the selected genes; a gene list; a parameter frame and a heatmap
#'
#' @import magrittr
#' @import dplyr
#' @import methods
#' @import viridis
#'
#' @export
#'
#' @examples  \dontrun{calMarker(expr = expr,thresh.1 = 0.5, thresh.2 = 0.4, num = 2, method = "del_MI")}
#'
searchMarker <- function(x,
                         thresh.1 = 0.5,
                         thresh.2 = NULL,
                         method = "del_MI",
                         num = 2,
                         gene.use = NULL,
                         meta.data = NULL,
                         ident.use = NULL,
                         scale.method = "scale",
                         lim.scale = 2,
                         border_color = "black",
                         colors = viridis(10)){
 UseMethod("searchMarker")
}

#' @import magrittr
#' @import dplyr
#' @import viridis
#' @export
searchMarker.matrix <- function(x,
                                thresh.1 = 0.5,
                                thresh.2 = NULL,
                                method = "del_MI",
                                num = 2,
                                gene.use = NULL,
                                meta.data = NULL,
                                ident.use = NULL,
                                scale.method = "scale",
                                lim.scale = 2,
                                border_color = "black",
                                colors = viridis(10)){
  #calculating average expression matrix
  if(ncol(x) >= 100){warning("more than 100 clusters detected... please check the input data")}

  expr.use <- x

  clusters <- colnames(expr.use) %>% as.vector()
  n.cluster <- length(clusters)


  #max_normalize the data
  message("max normalizing the matirx...")
  expr.norm <- apply(t(expr.use),2,function(x){x/max(x)}) %>% t()


  #maximumly expressed cluster
  expr.norm <- as.data.frame(expr.norm)
  expr.norm$max.X <- apply(expr.norm, 1, function(x){clusters[which.max(x)]})
  expr.norm$mean.orig <- apply(expr.use, 1, mean)


  #flitrate clusters accroding to thresh.2
  if(is.null(thresh.2)){
    message("setting thresh.2 as 0.1 by default")
    thresh.2 <- 0.1
  }


  message("screen genes in each cluster according to thresh.2")
  expr.norm <- filtGene(expr.norm,thresh.2 = thresh.2, clusters = clusters)


  # calculating parameters of mean.2 mean.3 NMI PMI...
  message("calculating parameters...")
  mean.frame <- ParaFrame(expr.norm,n.cluster,thresh.1,clusters)


  # calculate method
  mean.frame$del_MI <- mean.frame$PMI - mean.frame$NMI


  # range
  mean.frame <- mean.frame[order(mean.frame[,"max.X"],mean.frame[,"n"],-mean.frame[,method]),]
  genes.markers <- with(mean.frame,
                        by(mean.frame[,"gene"],
                           max.X,function(x){head(x,num)})) %>% unlist() %>% as.character() %>% rev()

  #plot heatmap
  message("Using \"RNA\" as the assay to plot Heatmap...")
  p <- HeatPlot(x = expr.use,
                genes = genes.markers,
                scale.method = scale.method,
                lim.scale = lim.scale,
                border_color = border_color,
                colors = colors)

  # a list containing all the results
  message("createing out put data...")
  calmarkers.out <- list(
    para_frame = mean.frame,
    genes.markers = genes.markers,
    exprs.markers = expr.use[genes.markers,],
    heatmap = p
  )
  message("Analysing Complete!")


  return(calmarkers.out)
  #warning 由于不是seurat对象，我们无法进行堆叠小提琴图的绘制
}

#' @import Seurat
#' @import magrittr
#' @import dplyr
#' @import viridis
#' @export
searchMarker.Seurat <- function(x,
                                thresh.1 = 0.5,
                                thresh.2 = NULL,
                                method = "del_MI",
                                num = 2,
                                gene.use = NULL,
                                meta.data = NULL,
                                ident.use = NULL,
                                scale.method = "scale",
                                lim.scale = 2,
                                border_color = "black",
                                colors = viridis(10)){

  #stop
  if(length(Idents(x)) == 1){stop("Ident with only one level, starTracer is not able to define the marker gene...")}

  #calculating average expression matrix
  .idents <- Seurat::Idents(x) %>% levels() %>% paste0(collapse = ", ")
  message(paste0("using ",.idents," as ident..."))

  .tmp <- any(grepl("FindVariableFeatures",Command(x)))

  if(is.null(gene.use)){
    message("using all genes as input features...")
    expr.use <- Seurat::AverageExpression(x)[[1]]
  } else if(gene.use == "HVG"){
    if(.tmp){
      message("using HVG as input features...")
      x <- x[VariableFeatures(x),]
      expr.use <- Seurat::AverageExpression(x)[[1]]
      } else {stop("please run Seurat::FindVariableFeatures()")}
  } else {
    message("plese set gene.use as \"HVG\" to use only HUV or NULL to use all genes")
  }

  clusters <- colnames(expr.use) %>% as.vector()
  n.cluster <- length(clusters)


  #max_normalize the data
  message("max normalizing the matirx...")
  expr.norm <- apply(t(expr.use),2,function(x){x/max(x)}) %>% t()


  #maximumly expressed cluster
  expr.norm <- as.data.frame(expr.norm)
  expr.norm$max.X <- apply(expr.norm, 1, function(x){clusters[which.max(x)]})
  expr.norm$mean.orig <- apply(expr.use, 1, mean)


  #flitrate clusters accroding to thresh.2
  if(is.null(thresh.2)){
    message("setting thresh.2 as 0.1 by default")
    thresh.2 <- 0.1
  }


  message("screen genes in each cluster according to thresh.2")
  expr.norm <- filtGene(expr.norm,thresh.2 = thresh.2, clusters = clusters)


  # calculating parameters of mean.2 mean.3 NMI PMI...
  message("calculating parameters...")
  mean.frame <- ParaFrame(expr.norm,n.cluster,thresh.1,clusters)


  # calculate method
  mean.frame$del_MI <- mean.frame$PMI - mean.frame$NMI


  # range
  mean.frame <- mean.frame[order(mean.frame[,"max.X"],mean.frame[,"n"],-mean.frame[,method]),]
  genes.markers <- with(mean.frame,
                        by(mean.frame[,"gene"],
                           max.X,function(x){head(x,num)})) %>% unlist() %>% as.character() %>% rev()


  #plot heatmap
  message("Using \"RNA\" as the assay to plot Heatmap...")
  p <- HeatPlot(x = expr.use,
                genes = genes.markers,
                scale.method = scale.method,
                lim.scale = lim.scale,
                border_color = border_color,
                colors = colors)


  # a list containing all the results
  message("createing out put data...")
  calmarkers.out <- list(
    para_frame = mean.frame,
    genes.markers = genes.markers,
    exprs.markers = expr.use[genes.markers,],
    heatmap = p
  )
  message("Analysing Complete!")
  return(calmarkers.out)
}

#' @import Matrix
#' @import magrittr
#' @import dplyr
#' @import pheatmap
#' @import stats
#' @import viridis
#' @export
searchMarker.dgCMatrix <- function(x,
                                   thresh.1 = 0.5,
                                   thresh.2 = NULL,
                                   method = "del_MI",
                                   num = 2,
                                   gene.use = NULL,
                                   meta.data = NULL,
                                   ident.use = NULL,
                                   scale.method = "scale",
                                   lim.scale = 2,
                                   border_color = "black",
                                   colors = viridis(10)){

  #get meta.data, store as "data":
  data <- data.frame(row.names = row.names(meta.data), idents = meta.data[,ident.use])
  if(class(meta.data[,ident.use]) == "factor"){
    data$idents <- factor(data$idents,levels = levels(meta.data[,ident.use]))
  } else {
    warning("metadata is not a factor, setting as one...")
    data$idents <- factor(data$idents)
  }

  #stop and message
  if(length(unique(data$idents)) == 1){stop("Ident with only one level, starTracer is not able to define the marker gene...")}
  message(paste0("using ",unique(data$idents)," as ident..."))

  #processing dgCMatrix

  category.matrix <- sparse.model.matrix(object = as.formula(
    object = paste0(
      '~0+',
      paste0(
        "data[,",
        1:length(x = ident.use),
        "]",
        collapse = ":"
        )
      )
    )
  )
  colsums <- colSums(x = category.matrix)
  # filter out cells without gene expression
  category.matrix <- category.matrix[, colsums > 0]
  colsums <- colsums[colsums > 0]
  # calculate rate
  category.matrix <- sweep(
    x = category.matrix,
    MARGIN = 2,
    STATS = colsums,
    FUN = "/")
  # using expm1 to deal with matrix
  data.use <- expm1(x = x)
  # getting average exprssion matrix
  expr.use <- as.matrix(x = (data.use %*% category.matrix))
  colnames(expr.use) <- levels(meta.data[,ident.use])

  if(is.null(gene.use)){
    message("using all genes as input features...")
    expr.use <- expr.use
  }

  clusters <- colnames(expr.use) %>% as.vector()
  n.cluster <- length(clusters)


  #max_normalize the data
  message("max normalizing the matirx...")
  expr.norm <- apply(t(expr.use),2,function(x){x/max(x)}) %>% t()


  #maximumly expressed cluster
  expr.norm <- as.data.frame(expr.norm)
  expr.norm$max.X <- apply(expr.norm, 1, function(x){clusters[which.max(x)]})
  expr.norm$mean.orig <- apply(expr.use, 1, mean)


  #flitrate clusters accroding to thresh.2
  if(is.null(thresh.2)){
    message("setting thresh.2 as 0.5 by default")
    thresh.2 <- 0.5
  }


  message("screen genes in each cluster according to thresh.2")
  expr.norm <- filtGene(expr.norm,thresh.2 = thresh.2, clusters = clusters)


  # calculating parameters of mean.2 mean.3 NMI PMI...
  message("calculating parameters...")
  mean.frame <- ParaFrame(expr.norm,n.cluster,thresh.1,clusters)


  # calculate method
  mean.frame$del_MI <- mean.frame$PMI - mean.frame$NMI


  # range
  mean.frame <- mean.frame[order(mean.frame[,"max.X"],mean.frame[,"n"],-mean.frame[,method]),]
  genes.markers <- with(mean.frame,
                        by(mean.frame[,"gene"],
                           max.X,function(x){head(x,num)})) %>% unlist() %>% as.character() %>% rev()


  #plot heatmap
  message("Using \"RNA\" as the assay to plot Heatmap...")
  p <- HeatPlot(x = expr.use,
                genes = genes.markers,
                scale.method = scale.method,
                lim.scale = lim.scale,
                border_color = border_color,
                colors = colors)

  # a list containing all the results
  message("createing out put data...")
  calmarkers.out <- list(
    para_frame = mean.frame,
    genes.markers = genes.markers,
    exprs.markers = expr.use[genes.markers,],
    heatmap = p
  )
  message("Analysing Complete!")
  return(calmarkers.out)
}







