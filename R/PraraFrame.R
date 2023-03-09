#' ParaFrame
#'
#' @param expr.norm max normalized function calculated before, but with a extra line of max.X
#' @param n.cluster cluster number
#' @param thresh.1 threshold to tell high expression data from the low expression level
#' @param clusters clusters from the single cell experiment
#'
#' @return a dataframe of all the parameters
#'
#' @examples \dontrun{ParaFrame(expr.norm = expr.norm,n.cluster = n.cluster,thresh.1 = thresh.1)}
ParaFrame <- function(expr.norm,n.cluster,thresh.1,clusters){
  #N and n
  n <- apply(expr.norm,1,function(x){
    line <- as.numeric(x[1:n.cluster])
    n <- length(line[line > thresh.1])
    return(n)
  })
  N <- n.cluster

  #mean3 and mean2
  mean.3 <- apply(expr.norm,1,function(x){
    line <- as.numeric(x[1:n.cluster])
    mean <- mean(line[line < thresh.1])
    return(mean)
  })
  mean.2 <- apply(expr.norm,1,function(x){
    line <- as.numeric(x[1:n.cluster])
    mean <- mean(line[line >= thresh.1])
    return(mean)
  })

  #calculate max.X again: after filtering with thresh.2, max.X need to be calculated again
  max.X <- apply(expr.norm[,1:n.cluster], 1,
                 function(x){colnames(expr.norm[,1:n.cluster])[which.max(x)]})

  #rest_cluster
  rest_cluster <- apply(expr.norm, 1, function(x){
    line <- as.numeric(x[1:n.cluster])
    colnames(expr.norm)[1:n.cluster][line > thresh.1]
  } %>% paste(collapse = "_"))

  #calculating marker.index
  PMI <- (1-mean.2)/(1-thresh.1)
  NMI <- mean.3/thresh.1

  # make a data frame
  mean.frame <- data.frame(mean.2 = mean.2,mean.3 = mean.3,PMI = PMI, NMI = NMI,
                           #opt.thresh.1 = opt.thresh.1,shift.thresh.1 = shift.thresh.1,
                           max.X = max.X,rest_cluster = rest_cluster,
                           #sum_non0 = sum_non0, sum_all, var_non0 = var_non0,
                           n = n,N = N)

  mean.frame$max.X <- factor(mean.frame$max.X,levels = clusters)
  mean.frame$gene <- rownames(mean.frame)

  return(mean.frame)
}
