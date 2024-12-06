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

#' Heat Kernel Function
#'
#' @description
#' Computes the heat kernel, which is a radial basis function commonly used in machine learning
#' and graph-based algorithms. The heat kernel measures similarity between points based on
#' their Euclidean distance, with a bandwidth parameter controlling the scale of similarity.
#'
#' @details
#' The heat kernel is defined as:
#' \deqn{K(x) = exp(-x^2/(2\sigma^2))}
#' where x is the distance between points and \eqn{\sigma} is the bandwidth parameter.
#' The kernel assigns higher values to points that are closer together and lower values
#' to points that are further apart, with the rate of decay controlled by sigma.
#'
#' @param x A numeric vector or matrix of distances between points
#' @param sigma Positive numeric value specifying the bandwidth parameter (default: 1).
#'   Larger values result in slower decay of similarity with distance.
#'
#' @return A numeric vector or matrix of the same dimensions as x containing the computed
#'   heat kernel values in the range [0,1].
#'
#' @examples
#' # Simple example with a sequence of distances
#' x <- seq(0, 5, length.out = 100)
#' y <- heat_kernel(x, sigma = 1)
#' plot(x, y, type = "l", main = "Heat Kernel", xlab = "Distance", ylab = "Similarity")
#'
#' # Example with a distance matrix
#' X <- matrix(rnorm(20), ncol = 2)
#' D <- as.matrix(dist(X))
#' K <- heat_kernel(D, sigma = 0.5)
#'
#' @seealso
#' \code{\link{estimate_sigma}} for estimating an appropriate bandwidth parameter
#' \code{\link{normalized_heat_kernel}} for a normalized version of the heat kernel
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


#' Graph-based Weight Matrix Construction
#'
#' @description
#' Constructs a weight matrix representing similarities between instances in a dataset
#' using various neighborhood and weighting schemes. This is particularly useful
#' for creating adjacency matrices for graph-based learning algorithms.
#'
#' @details
#' The function provides several methods for constructing the weight matrix:
#'
#' \strong{Neighbor Modes:}
#' \itemize{
#'   \item \code{knn}: k-nearest neighbors approach
#'   \item \code{epsilon}: epsilon-neighborhood approach (not yet implemented)
#' }
#'
#' \strong{Weight Modes:}
#' \itemize{
#'   \item \code{heat}: Heat kernel weights
#'   \item \code{normalized}: Normalized heat kernel weights
#'   \item \code{binary}: Binary weights (1 for neighbors, 0 otherwise)
#'   \item \code{euclidean}: Euclidean distances as weights
#'   \item \code{cosine}: Cosine similarity weights
#'   \item \code{correlation}: Correlation-based weights
#' }
#'
#' @param X Numeric matrix where rows are instances and columns are features
#' @param k Integer specifying the number of nearest neighbors (default: 5)
#' @param neighbor_mode Character string specifying the neighbor selection method:
#'   "knn" or "epsilon"
#' @param weight_mode Character string specifying the weight computation method:
#'   "heat", "normalized", "binary", "euclidean", "cosine", or "correlation"
#' @param type Character string specifying the neighbor relationship type:
#'   "normal" (standard k-NN), "mutual" (mutual k-NN), or "asym" (asymmetric)
#' @param sigma Numeric bandwidth parameter for heat kernel (if NULL, estimated from data)
#' @param eps Numeric epsilon value for epsilon-neighborhood (not implemented)
#' @param labels Optional factor vector of instance labels for supervised methods
#' @param ... Additional arguments passed to internal functions
#'
#' @return An object of class "neighbor_graph" containing:
#' \itemize{
#'   \item A sparse weight matrix representing the graph
#'   \item Parameters used to construct the graph
#' }
#'
#' @examples
#' # Generate example data
#' set.seed(123)
#' X <- matrix(rnorm(100 * 10), nrow = 100)
#'
#' # Create graph with heat kernel weights
#' g1 <- graph_weights(X, k = 5, neighbor_mode = "knn",
#'                    weight_mode = "heat")
#'
#' # Create graph with binary weights and mutual neighbors
#' g2 <- graph_weights(X, k = 5, neighbor_mode = "knn",
#'                    weight_mode = "binary", type = "mutual")
#'
#' # Create graph with normalized weights and labels
#' labels <- factor(rep(1:2, each = 50))
#' g3 <- graph_weights(X, k = 5, neighbor_mode = "knn",
#'                    weight_mode = "normalized", labels = labels)
#'
#' @seealso
#' \code{\link{heat_kernel}} for the heat kernel function
#' \code{\link{estimate_sigma}} for bandwidth parameter estimation
#'
#' @references
#' Luxburg, U. (2007). A tutorial on spectral clustering.
#' Statistics and Computing, 17(4), 395-416.
#'
#' @export
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
#' @param ef Search parameter that controls speed/accuracy trade-off (default: 200).
#' @param M Parameter that controls index size/speed trade-off (default: 16).
#'
#' @return If 'as' is "index_sim", a two-column matrix where the first column contains the indices of nearest neighbors and the second column contains the corresponding kernel values.
#'         If 'as' is "igraph", an igraph object representing the cross adjacency graph.
#'         If 'as' is "sparse", a sparse adjacency matrix.
#'
#' @importFrom RcppHNSW hnsw_knn
#' @export
cross_adjacency <- function(X, Y, k=5, FUN=heat_kernel, type=c("normal", "mutual", "asym"),
                          as=c("igraph", "sparse", "index_sim"), ef=200, M=16) {

  assert_that(k > 0 && k <= nrow(X))

  as <- match.arg(as)
  type <- match.arg(type)

  # Build index and search using RcppHNSW
  nn <- RcppHNSW::hnsw_knn(X, query = Y, k = k + 1, distance = "l2", ef = ef, M = M)
  
  # RcppHNSW returns squared L2 distances, so we need to take sqrt
  # Also add small epsilon to avoid numerical issues
  nnd <- sqrt(nn$dist[, 2:ncol(nn$dist), drop=FALSE] + 1e-16)
  nni <- nn$idx[, 2:ncol(nn$idx), drop=FALSE]
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

#' Compute weighted k-nearest neighbors using HNSW
#'
#' @param X A matrix where each row is a data point
#' @param k Number of nearest neighbors (default: 5)
#' @param FUN Function to compute weights from distances (default: heat_kernel)
#' @param type Type of neighborhood relationship: "normal", "mutual", or "asym"
#' @param as Output format: "igraph" or "sparse"
#' @param ... Additional parameters passed to FUN
#' @importFrom RcppHNSW hnsw_build hnsw_search
#' @importFrom Matrix sparseMatrix
#' @export
weighted_knn <- function(X, k=5, FUN=heat_kernel,
                        type=c("normal", "mutual", "asym"),
                        as=c("igraph", "sparse"),
                        ...) {
  
  assert_that(k > 0 && k <= nrow(X))
  type <- match.arg(type)
  as <- match.arg(as)
  
  # Build HNSW index
  ann <- hnsw_build(X, distance = "l2", M = 16, ef = 100)
  
  # Search for k nearest neighbors
  nn_results <- hnsw_search(X, ann, k = k + 1)  # +1 because first neighbor is self
  
  # Remove self-connections (first column)
  nn_index <- nn_results$idx[, -1, drop=FALSE]
  nn_dist <- sqrt(nn_results$dist[, -1, drop=FALSE])  # Convert from squared L2 to Euclidean
  
  # Apply weight function to distances
  hval <- apply(nn_dist, 2, FUN, ...)
  
  # Convert to sparse matrix format
  n <- nrow(X)
  M <- indices_to_sparse(nn_index, hval, return_triplet=TRUE, idim=n, jdim=n)
  
  if (type == "mutual") {
    # For mutual kNN, only keep edges that appear in both directions
    M_sparse <- triplet_to_matrix(M, c(n,n))
    M_t <- t(M_sparse)
    M_mutual <- M_sparse * (M_sparse > 0 & M_t > 0)
    M_mutual <- (M_mutual + t(M_mutual))/2
    
    if (as == "igraph") {
      M_mutual <- as(M_mutual, "dgTMatrix")
      ret <- graph_from_adjacency_matrix(M_mutual, mode="undirected", weighted=TRUE)
    } else {
      ret <- M_mutual
    }
  } else if (type == "normal") {
    # For normal kNN, make it symmetric by averaging with transpose
    M_sparse <- triplet_to_matrix(M, c(n,n))
    M_sym <- (M_sparse + t(M_sparse))/2
    
    if (as == "igraph") {
      M_sym <- as(M_sym, "dgTMatrix")
      ret <- graph_from_adjacency_matrix(M_sym, mode="undirected", weighted=TRUE)
    } else {
      ret <- M_sym
    }
  } else {
    # For asymmetric kNN, keep it as is
    M_sparse <- triplet_to_matrix(M, c(n,n))
    
    if (as == "igraph") {
      M_sparse <- as(M_sparse, "dgTMatrix")
      ret <- graph_from_adjacency_matrix(M_sparse, mode="directed", weighted=TRUE)
    } else {
      ret <- M_sparse
    }
  }
  
  ret
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



#' Estimate Optimal Bandwidth Parameter for Heat Kernel
#'
#' @description
#' Estimates an appropriate bandwidth parameter (sigma) for the heat kernel
#' based on the distribution of pairwise distances in the data. This is a crucial
#' parameter that affects the scale at which similarities are measured.
#'
#' @details
#' The function estimates sigma by:
#' 1. Computing pairwise distances between points (or a sample of points)
#' 2. Finding the specified quantile of non-zero distances
#' This approach helps choose a sigma value that's appropriate for the scale
#' of the data.
#'
#' @param X Numeric matrix where rows are instances and columns are features
#' @param prop Numeric value in (0,1) specifying the quantile to use (default: 0.25)
#' @param nsamples Integer specifying maximum number of samples to use (default: 500)
#' @param normalized Logical indicating whether distances should be computed on
#'   normalized data (default: FALSE)
#'
#' @return Numeric value representing the estimated optimal bandwidth parameter
#'
#' @examples
#' # Generate example data
#' set.seed(123)
#' X <- matrix(rnorm(1000 * 10), nrow = 1000)
#'
#' # Estimate sigma with default parameters
#' sigma1 <- estimate_sigma(X)
#'
#' # Estimate sigma with different quantile
#' sigma2 <- estimate_sigma(X, prop = 0.1)
#'
#' # Estimate sigma for normalized data
#' sigma3 <- estimate_sigma(X, normalized = TRUE)
#'
#' @export
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
