

#' @keywords internal
as_triplet <- function(M) {
  tm <- as(M, "dgTMatrix")
  cbind(i=tm@i+1, j=tm@j+1, x=tm@x)
}

#' @keywords internal
triplet_to_matrix <- function(trip, dim) {
  sparseMatrix(i=trip[,1], j=trip[,2], x=trip[,3],dims=dim)
}


#' @keywords internal
indices_to_sparse <- function(nn.index, hval, return_triplet=FALSE,
                              idim=nrow(nn.index),
                              jdim=nrow(nn.index)) {


  M <- do.call(rbind, lapply(1:nrow(nn.index), function(i) {
    cbind(i, nn.index[i,], hval[i,])
  }))

  if (return_triplet) {
    M
  } else {
    Matrix::sparseMatrix(i=M[,1], j=M[,2], x=M[,3], dims=c(idim, jdim))
  }
}

#' Compute the Heat Kernel
#'
#' This function computes the heat kernel, which is a radial basis function that can be used for smoothing, interpolation, and approximation tasks. The heat kernel is defined as exp(-x^2/(2*sigma^2)), where x is the distance and sigma is the bandwidth. It acts as a similarity measure for points in a space, assigning high values for close points and low values for distant points.
#'
#' @section Details:
#' The heat kernel is widely used in various applications, including machine learning, computer graphics, and image processing. It can be employed in kernel methods, such as kernel PCA, Gaussian process regression, and support vector machines, to capture the local structure of the data. The heat kernel's behavior is controlled by the bandwidth parameter sigma, which determines the smoothness of the resulting function.
#'
#' @param x A numeric vector or matrix representing the distances between data points.
#' @param sigma The bandwidth of the heat kernel, a positive scalar value. Default is 1.
#'
#' @return A numeric vector or matrix with the same dimensions as the input `x`, containing the computed heat kernel values.
#'
#' @examples
#' x <- seq(-3, 3, length.out = 100)
#' y <- heat_kernel(x, sigma = 1)
#' plot(x, y, type = "l", main = "Heat Kernel")
#'
#' @export
heat_kernel <- function(x, sigma=1) {
  exp((-x^2)/(2*sigma^2))
}
#' inverse_heat_kernel
#'
#' @param x the distances
#' @param sigma the bandwidth
#' @export
inverse_heat_kernel <- function(x, sigma=1) {
  #exp((-x^2)/(2*sigma^2))
  exp((-2*sigma^2)/x)
}




#' normalized_heat_kernel
#'
#' @param x the distances
#' @param sigma the bandwidth
#' @param len the normalization factor (e.g. the length of the feature vectors)
#' @export
normalized_heat_kernel <- function(x, sigma=.68, len) {
  norm_dist <- (x^2)/(2*len)
  exp(-norm_dist/(2*sigma^2))
}


#' @keywords internal
correlation_kernel <- function(x, len) {
  1 - (x^2)/(2*(len-1))
}


#' @keywords internal
cosine_kernel <- function(x, sigma=1) {
  1 - (x^2)/2
}



#' @keywords internal
get_neighbor_fun <- function(weight_mode=c("heat", "binary", "normalized", "euclidean", "cosine", "correlation"),
                             len,sigma) {

  weight_mode <- match.arg(weight_mode)

  wfun <- if (weight_mode == "heat") {
    function(x) heat_kernel(x, sigma)
  } else if (weight_mode == "binary") {
    function(x) (x>0) * 1
  } else if (weight_mode == "normalized") {
    function(x) normalized_heat_kernel(x,sigma=sigma,len=len)
  } else if (weight_mode == "euclidean") {
    function(x) x
  } else if (weight_mode == "cosine") {
    function(x) cosine_kernel(x)
  } else if (weight_mode == "correlation") {
    function(x) correlation_kernel(x,len)
  }
}


#' Compute Similarity Matrix for Factors in a Data Frame
#'
#' Calculate the similarity matrix for a set of factors in a data frame using various similarity methods.
#'
#' @param des A data frame containing factors for which the similarity matrix will be computed.
#' @param method A character vector specifying the method used for computing the similarity. The available methods are:
#'   \itemize{
#'     \item "Jaccard" - Jaccard similarity coefficient
#'     \item "Rogers" - Rogers and Tanimoto similarity coefficient
#'     \item "simple matching" - Simple matching coefficient
#'     \item "Dice" - Dice similarity coefficient
#'   }
#' @return A similarity matrix computed using the specified method for the factors in the data frame.
#' @details
#' The \code{factor_sim} function computes the similarity matrix for a set of factors in a data frame using the chosen method.
#' The function first converts the data frame into a model matrix, then calculates the similarity matrix using the \code{proxy::simil}
#' function from the \code{proxy} package.
#'
#' The function supports four similarity methods: Jaccard, Rogers, simple matching, and Dice. The choice of method depends on the
#' specific use case and the desired properties of the similarity measure.
#'
#' @export
#' @examples
#' # Sample data
#' des <- data.frame(
#'   var1 = factor(c("a", "b", "a", "b", "a")),
#'   var2 = factor(c("c", "c", "d", "d", "d"))
#' )
#'
#' # Compute similarity matrix using Jaccard method
#' sim_jaccard <- factor_sim(des, method = "Jaccard")
#'
#' # Compute similarity matrix using Dice method
#' sim_dice <- factor_sim(des, method = "Dice")
factor_sim <- function(des, method=c("Jaccard", "Rogers", "simple matching", "Dice")) {
  Fmat <- do.call(cbind, lapply(1:ncol(des), function(i) {
    nam <- colnames(des)[i]
    model.matrix(as.formula(paste("~ ", nam, " -1")), data=des)
  }))

  proxy::simil(Fmat, method=method)
}


#' Compute Weighted Similarity Matrix for Factors in a Data Frame
#'
#' Calculate the weighted similarity matrix for a set of factors in a data frame.
#'
#' @param des A data frame containing factors for which the weighted similarity matrix will be computed.
#' @param wts A numeric vector of weights corresponding to the factors in the data frame. The default is equal weights for all factors.
#' @return A weighted similarity matrix computed for the factors in the data frame.
#'
#' @export
#' @examples
#' # Sample data
#' des <- data.frame(
#'   var1 = factor(c("a", "b", "a", "b", "a")),
#'   var2 = factor(c("c", "c", "d", "d", "d"))
#' )
#'
#' # Compute similarity matrix with default equal weights
#' sim_default_weights <- weighted_factor_sim(des)
#'
#' # Compute similarity matrix with custom weights
#' sim_custom_weights <- weighted_factor_sim(des, wts = c(0.7, 0.3))
weighted_factor_sim <- function(des, wts=rep(1, ncol(des))/ncol(des)) {
  wts <- wts/sum(wts)
  Fmat <- lapply(1:ncol(des), function(i) {
    nam <- colnames(des)[i]
    labs <- des[[nam]]
    label_matrix(labs, labs) * wts[i]
  })

  Reduce("+", Fmat)
}


#' Estimate Bandwidth Parameter (Sigma) for the Heat Kernel
#'
#' Estimate a reasonable bandwidth parameter (sigma) for the heat kernel based on a data matrix and the specified quantile of the frequency distribution of distances.
#'
#' @param X A data matrix where samples are rows and features are columns.
#' @param prop A numeric value representing the quantile of the frequency distribution of distances used to determine the bandwidth parameter. Default is 0.25.
#' @param nsamples An integer representing the number of samples to draw from the data matrix. Default is 500.
#' @param normalized A logical value indicating whether to normalize the data. Default is FALSE.
#'
#' @return A numeric value representing the estimated bandwidth parameter (sigma) for the heat kernel.
#' @export
#'
#' @examples
#' # Sample data
#' X <- matrix(rnorm(1000), nrow=100, ncol=10)
#'
#' # Estimate sigma with default parameters
#' sigma_default <- estimate_sigma(X)
#'
#' # Estimate sigma with custom quantile and number of samples
#' sigma_custom <- estimate_sigma(X, prop=0.3, nsamples=300)
estimate_sigma <- function(X, prop=.25, nsamples=500, normalized=FALSE) {
  ## estimate sigma
  if (nrow(X) <= 500) {
    nsamples <- nrow(X)
    sam <- 1:nrow(X)
  } else {
    nsamples <- min(nsamples, nrow(X))
    sam <- sample(1:nrow(X), nsamples)
  }

  d <- dist(X[sam,])

  qs <- quantile(d[d!=0],seq(prop, 1,by=.1))
  if (all(is.na(qs))) {
    stop("could not estimate sigma, all quantiles are NA ...")
  } else {
    qs <- qs[!is.na(qs)]
    qs[1]
  }
}


#' Convert a Data Matrix to an Adjacency Graph
#'
#' Convert a data matrix with n instances and p features to an n-by-n adjacency graph using specified neighbor and weight modes.
#'
#' @details
#' This function converts a data matrix with n instances and p features into an adjacency graph. The adjacency graph is created
#' based on the specified neighbor and weight modes. The neighbor mode determines how neighbors are assigned weights, and the weight
#' mode defines the method used to compute weights.
#'
#' @param X A data matrix where each row represents an instance and each column represents a variable. Similarity is computed over instances.
#' @param k An integer representing the number of neighbors (ignored when neighbor_mode is not 'epsilon').
#' @param neighbor_mode A character string specifying the method for assigning weights to neighbors, either "supervised", "knn", "knearest_misses", or "epsilon".
#' @param weight_mode A character string specifying the weight mode: binary (1 if neighbor, 0 otherwise), 'heat', 'normalized', 'euclidean', 'cosine', or 'correlation'.
#' @param type A character string specifying the nearest neighbor policy, one of: normal, mutual, asym.
#' @param sigma A numeric parameter for the heat kernel (exp(-dist/(2*sigma^2))).
#' @param eps A numeric value representing the neighborhood radius when neighbor_mode is 'epsilon' (not implemented).
#' @param labels A factor vector representing the class of the categories when weight_mode is 'supervised' with nrow(labels) equal to nrow(X).
#' @param ... Additional parameters passed to the internal functions.
#'
#' @return An adjacency graph based on the specified neighbor and weight modes.
#' @export
#'
#' @examples
#' # Sample data
#' X <- matrix(rnorm(100*100), 100, 100)
#' sm <- graph_weights(X, neighbor_mode="knn", weight_mode="normalized", k=3)
#'
#' labels <- factor(rep(letters[1:4], 5))
#' sm3 <- graph_weights(X, neighbor_mode="knn", k=3, labels=labels, weight_mode="cosine")
#' sm4 <- graph_weights(X, neighbor_mode="knn", k=100, labels=labels, weight_mode="cosine")
graph_weights <- function(X, k=5, neighbor_mode=c("knn", "epsilon"),
                          weight_mode=c("heat", "normalized", "binary", "euclidean",
                                        "cosine", "correlation"),
                          type=c("normal", "mutual", "asym"),
                          sigma=NULL,eps=NULL, labels=NULL, ...) {

  neighbor_mode = match.arg(neighbor_mode)
  weight_mode = match.arg(weight_mode)
  type <- match.arg(type)

  if (weight_mode == "normalized" || weight_mode == "correlation") {
    X <- t(scale(t(X), center=TRUE, scale=TRUE))
  } else if (weight_mode == "cosine") {
    X <- t(apply(X, 1, function(x) x/sqrt(sum(x^2))))
  }

  if ((is.null(sigma)) && (weight_mode %in% c("heat", "normalized"))) {
    if (weight_mode == "heat") {
      sigma <- estimate_sigma(X)
    } else if (weight_mode == "normalized") {
      sigma <- estimate_sigma(X, normalized=TRUE)
    }
    message("sigma is ", sigma)
  }

  wfun <- get_neighbor_fun(weight_mode, len=ncol(X), sigma=sigma)

  if (neighbor_mode == "knn") {
    W <- weighted_knn(X, k, FUN=wfun, type=type,...)
  } else if (neighbor_mode == "epsilon") {
    stop("epsilon not implemented")
  }

  neighbor_graph(W, params=list(k=k, neighbor_mode=neighbor_mode,
                                     weight_mode=weight_mode,
                                     sigma=sigma,
                                     type=type,
                                     labels=labels))



}



#' Threshold Adjacency
#'
#' This function extracts the k-nearest neighbors from an existing adjacency matrix.
#' It returns a new adjacency matrix containing only the specified number of nearest neighbors.
#'
#' @param A An adjacency matrix representing the graph.
#' @param k An integer specifying the number of neighbors to consider (default: 5).
#' @param type A character string indicating the type of k-nearest neighbors graph to compute. One of "normal" or "mutual" (default: "normal").
#' @param ncores An integer specifying the number of cores to use for parallel computation (default: 1).
#'
#' @return A sparse adjacency matrix containing only the specified number of nearest neighbors.
#'
#' @examples
#' A <- matrix(runif(100), 10, 10)
#' A_thresholded <- threshold_adjacency(A, k = 5)
#'
#' @importFrom assertthat assert_that
#' @importFrom parallel mclapply
#' @importFrom Matrix sparseMatrix
#' @export
threshold_adjacency <- function(A, k=5, type=c("normal", "mutual"), ncores=1) {
  assert_that(k > 0 && k <= nrow(A))

  type <- match.arg(type)

  jind <- 1:nrow(A)

  A2 <- do.call(rbind, parallel::mclapply(1:nrow(A), function(i) {
    ord <- order(A[i,], decreasing=TRUE)
    cbind(i=i, j=jind[ord[1:k]],x=A[i, ord[1:k]])
  }, mc.cores=ncores))

  m <- sparseMatrix(i=A2[,1], j=A2[,2], x=A2[,3], dims=c(nrow(A), nrow(A)))

  m <- if (type == "normal") {
    psparse(m, pmax)
  } else {
    psparse(m, pmin)
  }
}



#' Cross Adjacency
#'
#' This function computes the cross adjacency matrix or graph between two sets of points
#' based on their k-nearest neighbors and a kernel function applied to their distances.
#'
#' @param X A matrix of size nXk, where n is the number of data points and k is the dimensionality of the feature space.
#' @param Y A matrix of size pXk, where p is the number of query points and k is the dimensionality of the feature space.
#' @param k An integer indicating the number of nearest neighbors to consider (default: 5).
#' @param FUN A kernel function to apply to the Euclidean distances between data points (default: heat_kernel).
#' @param type A character string indicating the type of adjacency to compute. One of "normal", "mutual", or "asym" (default: "normal").
#' @param as A character string indicating the format of the output. One of "igraph", "sparse", or "index_sim" (default: "igraph").
#'
#' @return If 'as' is "index_sim", a two-column matrix where the first column contains the indices of nearest neighbors and the second column contains the corresponding kernel values.
#'         If 'as' is "igraph", an igraph object representing the cross adjacency graph.
#'         If 'as' is "sparse", a sparse adjacency matrix.
#' @export
cross_adjacency <- function(X, Y, k=5, FUN=heat_kernel, type=c("normal", "mutual", "asym"),
                          as=c("igraph", "sparse", "index_sim")) {

  assert_that(k > 0 && k <= nrow(X))

  as <- match.arg(as)
  type <- match.arg(type)

  nn <- rflann::Neighbour(as.matrix(Y), X, k=k+1)
  nnd <- sqrt(nn$distances[, 2:ncol(nn$distances),drop=FALSE] + 1e-16)
  nni <- nn$indices[, 2:ncol(nn$indices),drop=FALSE]
  hval <- FUN(nnd)

  if (as == "index_sim") {
    return(cbind(as.vector(nni), as.vector(hval)))
  }


  W <- if (k == 1) {
    indices_to_sparse(as.matrix(nni), as.matrix(hval), idim=nrow(X), jdim=nrow(X))
  } else {
    indices_to_sparse(nni, hval, idim=nrow(X), jdim=nrow(X))
  }

  gg <- if (type == "normal") {
    igraph::graph_from_adjacency_matrix(W, weighted=TRUE, mode="max")
  } else if (type == "mutual") {
    igraph::graph_from_adjacency_matrix(W, weighted=TRUE, mode="min")
  } else if (type == "asym") {
    igraph::graph_from_adjacency_matrix(W, weighted=TRUE, mode="directed")
  }

  if (as == "sparse") {
    igraph::as_adjacency_matrix(gg, attr="weight")
  } else {
    gg
  }

}

#' Weighted k-Nearest Neighbors
#'
#' This function computes a weighted k-nearest neighbors graph or adjacency matrix from a data matrix.
#' The function takes into account the Euclidean distance between instances and applies a kernel function
#' to convert the distances into similarities.
#'
#' @param X A data matrix where rows are instances and columns are features.
#' @param k An integer specifying the number of nearest neighbors to consider (default: 5).
#' @param FUN A kernel function used to convert Euclidean distances into similarities (default: heat_kernel).
#' @param type A character string indicating the type of k-nearest neighbors graph to compute. One of "normal", "mutual", or "asym" (default: "normal").
#' @param as A character string specifying the format of the output. One of "igraph" or "sparse" (default: "igraph").
#' @param ... Additional arguments passed to the nearest neighbor search function (rflann::Neighbour).
#'
#' @return If 'as' is "igraph", an igraph object representing the weighted k-nearest neighbors graph.
#'         If 'as' is "sparse", a sparse adjacency matrix.
#'
#' @examples
#' X <- matrix(rnorm(10 * 10), 10, 10)
#' w <- weighted_knn(X, k = 5)
#'
#' @importFrom assertthat assert_that
#' @importFrom igraph graph_from_adjacency_matrix as_adjacency_matrix
#' @importFrom Matrix sparseMatrix
#' @importFrom rflann Neighbour
#' @export
weighted_knn <- function(X, k=5, FUN=heat_kernel,
                         type=c("normal", "mutual", "asym"),
                         as=c("igraph", "sparse"),
                         ...) {
  assert_that(k > 0 && k <= nrow(X))
  as <- match.arg(as)


  type <- match.arg(type)
  #nn <- FNN::get.knn(X, k=k)
  #nn <- nabor::knn(X, k=k)
  nn <- rflann::Neighbour(X, X,k=min(k+1, nrow(X)), ...)
  #nnd <- nn$nn.dist + 1e-16

  nnd <- sqrt(nn$distances[, 2:ncol(nn$distances),drop=FALSE] + 1e-16)
  nni <- nn$indices[, 2:ncol(nn$indices), drop=FALSE]
  hval <- FUN(nnd)

  #W <- indices_to_sparse(nn$nn.index, hval)
  W <- if (k == 1) {
    indices_to_sparse(as.matrix(nni), as.matrix(hval), idim=nrow(X), jdim=nrow(X))
  } else {
    indices_to_sparse(nni, hval, idim=nrow(X), jdim=nrow(X))
  }

  ##bW <- W != 0
  ##tmp <- as(x$G, "dgTMatrix")
  ## should return the raw distances?

  gg <- if (type == "normal") {
    igraph::graph_from_adjacency_matrix(W, weighted=TRUE, mode="max")
  } else if (type == "mutual") {
    igraph::graph_from_adjacency_matrix(W, weighted=TRUE, mode="min")
  } else if (type == "asym") {
    igraph::graph_from_adjacency_matrix(W, weighted=TRUE, mode="directed")
  }

  if (as == "sparse") {
    igraph::as_adjacency_matrix(gg, attr="weight")
  } else {
    gg
  }
}

#' Apply a Function to Non-Zero Elements in a Sparse Matrix
#'
#' This function applies a specified function (e.g., max) to each pair of non-zero elements in a sparse matrix.
#' It can return the result as a triplet representation or a sparse matrix.
#'
#' @param M A sparse matrix object from the Matrix package.
#' @param FUN A function to apply to each pair of non-zero elements in the sparse matrix.
#' @param return_triplet A logical value indicating whether to return the result as a triplet representation. Default is FALSE.
#'
#' @return If return_triplet is TRUE, a matrix containing the i, j, and x values in the triplet format; otherwise, a sparse matrix with the updated values.
#'
#' @importFrom Matrix which
#' @importFrom Matrix sparseMatrix
#' @examples
#' library(Matrix)
#' M <- sparseMatrix(i = c(1, 3, 1), j = c(2, 3, 3), x = c(1, 2, 3))
#' psparse_max <- psparse(M, FUN = max)
#' psparse_sum_triplet <- psparse(M, FUN = `+`, return_triplet = TRUE)
#' @export
psparse <- function(M, FUN, return_triplet=FALSE) {
  ind <- which(M != 0, arr.ind=TRUE)
  x1 <- M[ind]
  x2 <- M[cbind(ind[,2], ind[,1])]

  if (return_triplet) {
    cbind(i=c(ind[,1],ind[,2]), j=c(ind[,2],ind[,1]), x=rep(FUN(x1,x2),2))
  } else {
    sm <- sparseMatrix(i=c(ind[,1],ind[,2]), j=c(ind[,2],ind[,1]), x=rep(FUN(x1,x2),2),
                       dims=dim(M), use.last.ij=TRUE)
    sm
  }
}



# psparse <- function(..., fun=c("max", "min"), na.rm=FALSE, return_triplet=FALSE) {
#   fun <- match.arg(fun)
#   # check that all matrices have conforming sizes
#   num.rows <- unique(sapply(list(...), nrow))
#   num.cols <- unique(sapply(list(...), ncol))
#   stopifnot(length(num.rows) == 1)
#   stopifnot(length(num.cols) == 1)
#
#   cat.summary <- do.call(rbind, lapply(list(...), summary))
#
#
#   out.summary <- if (fun == "min") {
#     aggregate(x ~ i + j, data = cat.summary, FUN=function(x) {
#       if (length(x) == 1) 0 else x[1]
#     })
#   } else {
#     aggregate(x ~ i + j, data = cat.summary, FUN=max, na.rm)
#   }
#
#   if (return_triplet) {
#     cbind(i=out.summary$i, j=out.summary$j, x=out.summary$x)
#   } else {
#     sparseMatrix(i = out.summary$i,
#                  j = out.summary$j,
#                  x = out.summary$x,
#                  dims = c(num.rows, num.cols))
#   }
# }
#
# pmin.sparse <- function(..., na.rm = FALSE, return_triplet=FALSE) { psparse(..., fun="min", return_triplet=return_triplet) }
#
# pmax.sparse <- function(..., na.rm = FALSE, return_triplet=FALSE) { psparse(..., fun="max", return_triplet=return_triplet ) }
#

