#' searchMarker
#'
#' this is the main function of searching marker genes
#'
#' @param x a Seurat Object, an average expression matrix or a sparse matrix.
#' @param thresh.1 the threshold for dividing clusters max-normalized gene expression to high or low, set to 0.5 by default
#' @param thresh.2 the the lowest expression level of the genes user prefers
#' @param thresh.2.neg the the highest expression level of the genes user prefers when calculating negative markers
#' @param method the arranging method of ordering the genes. Set default and strongly suggestted as "pos", starTracer will only find positivie marker genes, setting as "all" to out put both up-regulated and down-regulated markers
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
#' @examples  \dontrun{calMarker(expr = expr,thresh.1 = 0.5, thresh.2 = 0.4, num = 2, method = "pos")}
#'
searchMarker <- function(x,
                         thresh.1 = 0.5,
                         thresh.2 = NULL,
                         thresh.2.neg = NULL,
                         method = "pos",
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
                                thresh.2.neg = NULL,
                                method = "pos",
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

  #calculating negative marker gene
  if(method == "all") expr.norm.neg <- 1-expr.norm

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
  mean.frame <- mean.frame[order(mean.frame[,"max.X"],mean.frame[,"n"],-mean.frame[,"del_MI"]),]
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
    heatmap = p,
    expr.use = expr.use
  )
  message("Analysing Positive Marker Gene Complete!")

  #calculating negative markers
  if(method == "all") {
    #maximumly expressed cluster
    expr.norm.neg <- as.data.frame(expr.norm.neg)
    expr.norm.neg$max.X <- apply(expr.norm.neg, 1, function(x){clusters[which.max(x)]})
    expr.norm.neg$mean.orig <- apply(expr.use, 1, mean)


    #flitrate clusters accroding to thresh.2.neg
    if(is.null(thresh.2.neg)){
      message("setting thresh.2.neg as 0.5 by default")
      thresh.2.neg <- 0.5
    }


    message("screen genes in each cluster according to thresh.2.neg")
    expr.norm.neg <- filtGene.Neg(expr.norm.neg,thresh.2.neg = thresh.2.neg, clusters = clusters)


    # calculating parameters of mean.2 mean.3 NMI PMI...
    message("calculating parameters...")
    mean.frame <- ParaFrame(expr.norm.neg,n.cluster,thresh.1,clusters)


    # calculate method
    mean.frame$del_MI <- mean.frame$PMI - mean.frame$NMI


    # range
    mean.frame <- mean.frame[order(mean.frame[,"max.X"],mean.frame[,"n"],-mean.frame[,"PMI"]),]
    mean.frame <- mean.frame[!is.nan(mean.frame$PMI),]
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
    calmarkers.out.neg <- list(
      genes.markers.neg = genes.markers,
      exprs.markers.neg = expr.use[genes.markers.neg,],
      heatmap.neg = p
    )
    message("Analysing Negative Marker Gene Complete!")

    calmarkers.out <- append(calmarkers.out,calmarkers.out.neg)
  }

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
                                thresh.2.neg = NULL,
                                method = "pos",
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
    expr.use <- average.express(x, assay = "RNA", slot = "data", verbose = TRUE)
  } else if(gene.use == "HVG"){
    if(.tmp){
      message("using HVG as input features...")
      expr.use <- average.express(x, assay = "RNA", slot = "data", verbose = TRUE)
      expr.use <- expr.use[VariableFeatures(x),]
      } else {stop("please run Seurat::FindVariableFeatures()")}
  } else {
    message("plese set gene.use as \"HVG\" to use only HUV or NULL to use all genes")
  }

  clusters <- colnames(expr.use) %>% as.vector()
  n.cluster <- length(clusters)


  #max_normalize the data
  message("max normalizing the matirx...")
  expr.norm <- apply(t(expr.use),2,function(x){x/max(x)}) %>% t()

  if(method == "all") expr.norm.neg <- 1-expr.norm

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
  mean.frame <- mean.frame[order(mean.frame[,"max.X"],mean.frame[,"n"],-mean.frame[,"del_MI"]),]
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
    heatmap = p,
    expr.use = expr.use
  )
  message("Analysing Positive Marker Gene Completed!")

  #calculating negative markers
  if(method == "all") {
    #maximumly expressed cluster
    expr.norm.neg <- as.data.frame(expr.norm.neg)
    expr.norm.neg$max.X <- apply(expr.norm.neg, 1, function(x){clusters[which.max(x)]})
    expr.norm.neg$mean.orig <- apply(expr.use, 1, mean)


    #flitrate clusters accroding to thresh.2.neg
    if(is.null(thresh.2.neg)){
      message("setting thresh.2.neg as 0.5 by default")
      thresh.2.neg <- 0.5
    }


    message("screen genes in each cluster according to thresh.2.neg")
    expr.norm.neg <- filtGene.Neg(expr.norm.neg,thresh.2.neg = thresh.2.neg, clusters = clusters)


    # calculating parameters of mean.2 mean.3 NMI PMI...
    message("calculating parameters...")
    mean.frame <- ParaFrame(expr.norm.neg,n.cluster,thresh.1,clusters)


    # calculate method
    mean.frame$del_MI <- mean.frame$PMI - mean.frame$NMI


    # range
    mean.frame <- mean.frame[order(mean.frame[,"max.X"],mean.frame[,"n"],-mean.frame[,"PMI"]),]
    mean.frame <- mean.frame[!is.nan(mean.frame$PMI),]
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
    calmarkers.out.neg <- list(
      genes.markers.neg = genes.markers,
      exprs.markers.neg = expr.use[genes.markers,],
      heatmap.neg = p
    )
    message("Analysing Negative Marker Gene Completed!")

    calmarkers.out <- append(calmarkers.out,calmarkers.out.neg)
  }

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
                                   thresh.2.neg = NULL,
                                   method = "pos",
                                   num = 2,
                                   gene.use = NULL,
                                   meta.data = NULL,
                                   ident.use = NULL,
                                   scale.method = "scale",
                                   lim.scale = 2,
                                   border_color = "black",
                                   colors = viridis(10)){

  #stop and message
  if(length(unique(data$idents)) == 1){stop("Ident with only one level, starTracer is not able to define the marker gene...")}
  message(paste0("using ",unique(data$idents)," as ident..."))

  expr.use <- average.express.dgC(x,meta.data = meta.data,ident.use = ident.use)

  expr.norm <- apply(t(expr.use),2,function(x){x/max(x)}) %>% t()

  if(method == "all") expr.norm.neg <- 1-expr.norm
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
  mean.frame <- mean.frame[order(mean.frame[,"max.X"],mean.frame[,"n"],-mean.frame[,"del_MI"]),]
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
    heatmap = p,
    expr.use = expr.use
  )
  message("Analysing Positive Marker Gene Completed!")

  if(method == "all") {
    #maximumly expressed cluster
    expr.norm.neg <- as.data.frame(expr.norm.neg)
    expr.norm.neg$max.X <- apply(expr.norm.neg, 1, function(x){clusters[which.max(x)]})
    expr.norm.neg$mean.orig <- apply(expr.use, 1, mean)


    #flitrate clusters accroding to thresh.2.neg
    if(is.null(thresh.2.neg)){
      message("setting thresh.2.neg as 0.5 by default")
      thresh.2.neg <- 0.5
    }


    message("screen genes in each cluster according to thresh.2.neg")
    expr.norm.neg <- filtGene.Neg(expr.norm.neg,thresh.2.neg = thresh.2.neg, clusters = clusters)


    # calculating parameters of mean.2 mean.3 NMI PMI...
    message("calculating parameters...")
    mean.frame <- ParaFrame(expr.norm.neg,n.cluster,thresh.1,clusters)


    # calculate method
    mean.frame$del_MI <- mean.frame$PMI - mean.frame$NMI


    # range
    mean.frame <- mean.frame[order(mean.frame[,"max.X"],mean.frame[,"n"],-mean.frame[,"PMI"]),]
    mean.frame <- mean.frame[!is.nan(mean.frame$PMI),]
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
    calmarkers.out.neg <- list(
      genes.markers.neg = genes.markers,
      exprs.markers.neg = expr.use[genes.markers.neg,],
      heatmap.neg = p
    )
    message("Analysing Negative Marker Gene Completed!")

    calmarkers.out <- append(calmarkers.out,calmarkers.out.neg)
  }

  return(calmarkers.out)
}







