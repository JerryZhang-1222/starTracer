#' para_frame example data
#'
#' This is an example data set calculated by starTracer::searchMarker().
#' Here we present mean.2, mean.3, NMI, PMI and del_MI as the parameters for each gene
#'
#' @format A matrix with 10 variables:
#' \describe{
#' \item{mean.2}{the average max.normalized expression of clusters above thresh.1}
#' \item{mean.3}{the average max.normalized expression of clusters below thresh.1}
#' \item{PMI}{the positive index calculated from mean.2}
#' \item{NMI}{the negative index calculated from mean.2}
#' \item{max.X}{the cluster with the maximum expression of each gene}
#' \item{rest_cluster}{other clusters with the expressinon of such gene}
#' \item{n}{the number of clusters with the average max.nomalized expression above thresh.1}
#' \item{N}{the number of all clusters}
#' \item{gene}{SYMBOL of the gene}
#' \item{del_MI}{the delta MI}
#' }
"data_para_frame"
