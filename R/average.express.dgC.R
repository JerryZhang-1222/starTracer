#' Title This funciton calculate the average expression matrix from the sparse matrix.
#'
#' @param object a Seurat object
#' @param ident the column that will be used as identities from the meta data
#' @param assay the assay that will be used
#' @param slot the slot that will be used
#' @param verbose
#'
#' @return An average expression matrix
#'
#' @import Matrix
#' @import Seurat
#'
#' @examples \dontrun{average.express(x = dgCMatrix,meta.data = meta.data,ident.use = "avtive.ident", verbose = TRUE)}
average.express.dgC <- function(mat, meta.data,ident.use, verbose = TRUE) {

  # convert data to a sparse matrix if it isn't already
  ifelse (class(data)[1] != "dgCMatrix",data <- Matrix(mat, sparse = TRUE),
          data <- mat)

  # get the idents class over which to average
  idents <- meta.data[,ident.use]

  # loop through all idents, averaging them in data
  ident.names <- unique(idents)

  if (verbose > 0) pb <- txtProgressBar(char = "=", style = 3, max = length(ident.names), width = 50)
  m <- list()
  data <- expm1(x = data)
  for (i in 1:length(ident.names)) {
    ifelse(idents %in% ident.names[i] %>% base::sum() > 1,
           m[[i]] <- Matrix::rowMeans(data[, which(idents == ident.names[i])]),
           m[[i]] <- mean(data[, which(idents == ident.names[i])]))

    if (verbose > 0) setTxtProgressBar(pb = pb, value = i)
  }
  result <- do.call(cbind, m)
  colnames(result) <- ident.names
  return(result)
}
