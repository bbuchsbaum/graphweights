#' Create a Nearest Neighbor Searcher Object
#'
#' @description
#' Creates an efficient nearest neighbor search object using the Hierarchical Navigable Small World
#' (HNSW) algorithm. This implementation is particularly useful for high-dimensional data and
#' large-scale nearest neighbor search tasks.
#'
#' @details
#' The HNSW algorithm builds a multilayer graph structure for efficient approximate nearest
#' neighbor search. The algorithm provides a good balance between search accuracy and speed.
#' Key parameters that affect the search quality and speed are:
#' * M: The number of connections per element in the graph
#' * ef: The size of the dynamic candidate list during search
#'
#' The supported distance metrics are:
#' * l2: Euclidean distance (L2 norm)
#' * euclidean: Same as l2
#' * cosine: Cosine similarity
#' * ip: Inner product
#'
#' @param X Numeric matrix where rows represent observations and columns represent features
#' @param labels Vector of labels for the observations (default: 1:nrow(labels))
#' @param ... Additional arguments passed to the HNSW builder
#' @param distance Character string specifying the distance metric (default: "l2")
#' @param M Integer controlling the number of connections (default: 16)
#' @param ef Integer controlling the search accuracy (default: 200)
#'
#' @return An object of class "nnsearcher" containing:
#' * X: The input feature matrix
#' * labels: Vector of observation labels
#' * ann: The HNSW index object
#' * distance: The chosen distance metric
#' * ef: The ef parameter value
#' * M: The M parameter value
#'
#' @examples
#' # Create searcher for iris dataset
#' data(iris)
#' X <- as.matrix(iris[, 1:4])
#' labels <- iris$Species
#' searcher <- nnsearcher(X, labels, distance = "l2", M = 16, ef = 200)
#'
#' # Search for nearest neighbors
#' nn <- find_nn(searcher, k = 5)
#'
#' # Create similarity graph
#' graph <- neighbor_graph(searcher, k = 5, type = "normal", transform = "heat")
#'
#' @seealso
#' * \code{\link{find_nn}} for searching nearest neighbors
#' * \code{\link{neighbor_graph}} for creating similarity graphs
#' * \code{\link{dist_to_sim}} for converting distances to similarities
#'
#' @import RcppHNSW
#' @export
nnsearcher <- function(X, labels=1:nrow(labels), ...,
                       distance=c("l2", "euclidean", "cosine", "ip"), M=16, ef=200) {
  distance <- match.arg(distance)
  X <- as.matrix(X)
  ann <- RcppHNSW::hnsw_build(X, distance, ef=ef, M=M)

  structure(list(
    X=X,
    labels=labels,
    ann=ann,
    distance=distance,
    ef=ef,
    M=M),
    class="nnsearcher")
}

#' Format Search Results
#'
#' @description
#' Formats the results from a nearest neighbor search, adding metadata about the
#' feature dimensionality and distance metric used.
#'
#' @param x A nnsearcher object
#' @param result Raw search results from HNSW
#'
#' @return A formatted search result object of class "nn_search"
#'
#' @export
search_result.nnsearcher <- function(x, result) {
  attr(result, "len") <- ncol(x$X)
  attr(result, "metric") <- x$distance
  class(result) <- "nn_search"
  result
}

#' Convert Distances to Similarities
#'
#' @description
#' Converts distance measures to similarity scores using various transformation methods.
#'
#' @details
#' Available transformation methods:
#' * heat: exp(-d^2/sigma)
#' * binary: 1 if connected, 0 otherwise
#' * normalized: 1 - d/max(d)
#' * cosine: cos(d)
#' * correlation: correlation coefficient
#'
#' @param x An nn_search object containing distances
#' @param method Character string specifying the transformation method
#' @param sigma Scaling parameter for heat kernel transformation
#'
#' @return Modified nn_search object with distances converted to similarities
#'
#' @export
dist_to_sim.nn_search <- function(x, method = c("heat", "binary", "normalized", "cosine", "correlation"), sigma=1) {
  method <- match.arg(method)
  fun <- get_neighbor_fun(method, x$len, x$sigma)
  x$dist <- fun(x$dist)
  x
}

#' Convert Distance Matrix to Similarity Matrix
#'
#' @description
#' Converts a sparse distance matrix to a similarity matrix using various transformation methods.
#'
#' @param x A sparse Matrix object containing distances
#' @param method Character string specifying the transformation method
#' @param sigma Scaling parameter for heat kernel transformation
#' @param len Feature dimensionality for normalization
#'
#' @return Modified Matrix with distances converted to similarities
#'
#' @export
dist_to_sim.Matrix <- function(x, method = c("heat", "binary", "normalized", "cosine", "correlation"), sigma=1, len=1) {
  method <- match.arg(method)
  fun <- get_neighbor_fun(method, len, sigma)

  wh <- which(x != 0)
  v <- fun(x[wh])
  x[wh] <- v
  x
}

#' Convert Search Results to Adjacency Matrix
#'
#' @description
#' Creates a sparse adjacency matrix from nearest neighbor search results.
#'
#' @param x An nn_search object
#' @param idim Number of rows in output matrix
#' @param jdim Number of columns in output matrix
#' @param return_triplet Logical; if TRUE, returns triplet format
#'
#' @return A sparse adjacency matrix or triplet representation
#'
#' @export
adjacency.nnsearch <- function(x, idim=nrow(x$idx), jdim=max(x$idx), return_triplet=FALSE) {
  indices_to_sparse(as.matrix(x$idx), as.matrix(x$dist), return_triplet=FALSE, idim=idim, jdim=jdim)
}

#' Find Nearest Neighbors
#'
#' @description
#' Performs nearest neighbor search using a pre-built HNSW index.
#'
#' @param x A nnsearcher object
#' @param query Optional query matrix for searching external points
#' @param k Number of nearest neighbors to find
#'
#' @return An nn_search object containing:
#' * idx: Matrix of indices to nearest neighbors
#' * dist: Matrix of distances to nearest neighbors
#' * labels: Matrix of labels of nearest neighbors
#'
#' @examples
#' # Search in iris dataset
#' data(iris)
#' X <- as.matrix(iris[, 1:4])
#' searcher <- nnsearcher(X)
#' 
#' # Find 5 nearest neighbors for all points
#' nn <- find_nn(searcher, k = 5)
#' 
#' # Search for specific query points
#' query <- matrix(rnorm(8), 2, 4)
#' nn_query <- find_nn(searcher, query = query, k = 3)
#'
#' @export
find_nn.nnsearcher <- function(x, query=NULL, k=5) {
  ret <- if (!is.null(query)) {
    chk::chk_matrix(query)
    RcppHNSW::hnsw_search(query, x$ann, k = k)
  } else {
    RcppHNSW::hnsw_search(x$X, x$ann, k = k)
  }

  ret$labels <- t(apply(ret$idx, 1,  function(i) x$labels[i]))
  search_result(x,ret)
}

#' Find Nearest Neighbors Among Subset
#'
#' @description
#' Performs nearest neighbor search within a specified subset of points.
#'
#' @param x A nn_searcher object
#' @param k Number of nearest neighbors to find
#' @param idx Indices of points to search among
#'
#' @return An nn_search object for the subset
#'
#' @export
find_nn_among.nn_searcher <- function(x, k=5, idx) {
  chk::chk_numeric(idx)

  X1 <- x$X[idx,]
  ann <- RcppHNSW::hnsw_build(X1, x$distance, M=x$M, ef=x$ef)
  nnres <- RcppHNSW::hnsw_search(X1, ann, k=k)
  class(nnres) <- "nn_search"
  nnres
}

#' Find Nearest Neighbors Between Groups
#'
#' @description
#' Performs nearest neighbor search between two specified groups of points.
#'
#' @param x A nn_searcher object
#' @param k Number of nearest neighbors to find
#' @param idx1 Indices of first group
#' @param idx2 Indices of second group
#' @param restricted Logical; if TRUE, restricts search to specified indices
#'
#' @return An nn_search object containing between-group neighbors
#'
#' @export
find_nn_between.nn_searcher <- function(x, k=5, idx1, idx2, restricted=FALSE) {
  chk::chk_numeric(idx1)
  chk::chk_numeric(idx2)

  ret <- if (!restricted) {
    X1 <- x$X[idx1,]
    X2 <- x$X[idx2,]
    ann <- RcppHNSW::hnsw_build(X1, x$distance, M=x$M, ef=x$ef)
    nnres <- RcppHNSW::hnsw_search(X2, ann, k=k)
    nnres$idx <- do.call(rbind, lapply(1:nrow(nnres$idx), function(i) {
      idx1[nnres$idx[i,]]
    }))
    nnres$labels <- t(apply(nnres$idx, 1,  function(i) x$labels[i]))
    nnres
  } else {
    find_nn.nnsearcher(x, x$X[idx2,])
  }

  class(nnres) <- "nn_search"
  nnres
}

#' Create Neighbor Graph from Nearest Neighbor Search
#'
#' @description
#' Creates a neighbor graph from nearest neighbor search results with various
#' graph types and similarity transformations.
#'
#' @details
#' Available graph types:
#' * normal: Regular directed graph
#' * asym: Asymmetric directed graph
#' * mutual: Mutual nearest neighbor graph
#'
#' Similarity transformations:
#' * heat: Heat kernel
#' * binary: Binary weights
#' * euclidean: Raw distances
#' * normalized: Normalized distances
#' * cosine: Cosine similarity
#' * correlation: Correlation coefficient
#'
#' @param x A nn_searcher object
#' @param query Optional query points
#' @param k Number of neighbors per node
#' @param type Character string specifying graph type
#' @param transform Character string specifying similarity transformation
#' @param sigma Scaling parameter for heat kernel
#'
#' @return A neighbor_graph object
#'
#' @examples
#' # Create neighbor graph from iris dataset
#' data(iris)
#' X <- as.matrix(iris[, 1:4])
#' searcher <- nnsearcher(X)
#' 
#' # Create different types of graphs
#' g1 <- neighbor_graph(searcher, k = 5, type = "normal", transform = "heat")
#' g2 <- neighbor_graph(searcher, k = 5, type = "mutual", transform = "binary")
#'
#' @export
neighbor_graph.nn_searcher <- function(x, query=NULL, k=5, type=c("normal", "asym", "mutual"),
                                     transform=c("heat", "binary", "euclidean",
                                               "normalized", "cosine", "correlation"), sigma=1) {

  nn <- find_nn.nnsearcher(x,k=(k+1))
  D <- nn$dist[,2:(k+1)]

  nnd <- sqrt(nn$dist[, 2:ncol(nn$dist),drop=FALSE])
  nni <- nn$idx[, 2:ncol(nn$idx), drop=FALSE]

  hfun <- get_neighbor_fun(transform, len=ncol(x$X), sigma=sigma)
  hval <- hfun(D)

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

  neighbor_graph(gg, params=list(k=k,
                                transform=transform,
                                sigma=sigma,
                                type=type,
                                labels=x$labels))
}

#' @export 
find_nn_among.class_graph <- function(x, X, k=5) {
  searcher <- nnsearcher(X,x$labels)
  ret <- lapply(x$class_indices, function(ind) {
    nnr <- find_nn_among.nn_searcher(searcher, k=k, ind)
    nnr$idx <- do.call(rbind, lapply(1:nrow(nnr$idx), function(i) {
      ind[nnr$idx[i,]]
    }))

    nnr$labels <- t(apply(nnr$idx, 1,  function(i) x$labels[i]))
    nnr
  })

  idx <- do.call(rbind, lapply(ret, "[[", "idx"))
  d <- do.call(rbind, lapply(ret, "[[", "dist"))
  labels <- do.call(rbind, lapply(ret, "[[", "labels"))
  search_result(searcher, list(idx=idx, dist=d, labels=labels))
}
