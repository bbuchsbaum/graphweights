#' Nearest Neighbor Searcher
#'
#' Create a nearest neighbor searcher object using HNSW (Hierarchical Navigable 
#' Small World) algorithm for efficient approximate nearest neighbor search.
#'
#' @param X A numeric matrix where each row represents a data point.
#' @param labels A vector of labels corresponding to each row in X. Defaults to row indices.
#' @param ... Additional arguments (currently unused).
#' @param distance The distance metric to use. One of "l2", "euclidean", "cosine", or "ip".
#' @param M The maximum number of connections for every new element during construction.
#' @param ef The size of the dynamic candidate list.
#'
#' @return An object of class "nnsearcher" containing the data matrix, labels, 
#'   HNSW index, and search parameters.
#'
#' @examples
#' \donttest{
#' # Create example data
#' X <- matrix(rnorm(100), nrow=10, ncol=10)
#' searcher <- nnsearcher(X)
#' }
#'
#' @importFrom RcppHNSW hnsw_build hnsw_search
#' @importFrom chk chk_matrix chk_numeric
#' @export
nnsearcher <- function(X, labels=1:nrow(X), ...,
                       distance=c("l2", "euclidean", "cosine", "ip"), M=16, ef=200) {
  distance <- match.arg(distance)
  X <- as.matrix(X)
  
  # Validation checks
  stopifnot(
    "Number of labels must equal the number of rows in the data matrix 'X'" = 
      length(labels) == nrow(X)
  )
  if (anyNA(labels)) {
    stop("NA values are not permitted in 'labels'.", call. = FALSE)
  }
  
  ann <- hnsw_build(X, distance, ef=ef, M=M)

  structure(list(
    X=X,
    labels=labels,
    ann=ann,
    distance=distance,
    ef=ef,
    M=M),
    class="nnsearcher")
}

#' Convert Search Results for nnsearcher Objects
#'
#' Convert raw search results to a standardized format for nnsearcher objects.
#'
#' @param x An object of class "nnsearcher".
#' @param result A raw result object from nearest neighbor search.
#'
#' @return An object of class "nn_search" with standardized field names.
#'
#' @method search_result nnsearcher
#' @export
search_result.nnsearcher <- function(x, result) {
  # Map field names to public API standard
  result$indices <- result$idx
  result$distances <- result$dist
  result$idx <- NULL      # Remove internal names
  result$dist <- NULL
  
  attr(result, "len") <- ncol(x$X)
  attr(result, "metric") <- x$distance
  class(result) <- "nn_search"
  result
}

#' Convert Distance to Similarity for nn_search Objects
#'
#' Convert distance values in a nearest neighbor search result to similarity values.
#'
#' @param x An object of class "nn_search".
#' @param method The transformation method for converting distances to similarities.
#' @param sigma The bandwidth parameter for the heat kernel method.
#'
#' @return The modified nn_search object with distances converted to similarities.
#'
#' @method dist_to_sim nn_search
#' @export
dist_to_sim.nn_search <- function(x, method = c("heat", "binary", "normalized", "cosine", "correlation"), sigma=1) {
  method <- match.arg(method)
  fun <- get_neighbor_fun(method, x$len, x$sigma)
  x$dist <- fun(x$dist)
  x
}


#' Convert Distance to Similarity for Matrix Objects
#'
#' Convert distance values in a sparse Matrix to similarity values.
#'
#' @param x A Matrix object containing distances.
#' @param method The transformation method for converting distances to similarities.
#' @param sigma The bandwidth parameter for the heat kernel method.
#' @param len The length parameter used in transformation calculations.
#'
#' @return The Matrix object with distances converted to similarities.
#'
#' @method dist_to_sim Matrix
#' @export
dist_to_sim.Matrix <- function(x, method = c("heat", "binary", "normalized", "cosine", "correlation"), sigma=1, len=1) {
  method <- match.arg(method)
  fun <- get_neighbor_fun(method, len, sigma)

  wh <- which(x != 0)
  v <- fun(x[wh])
  x[wh] <- v
  x
}

#' Create Adjacency Matrix from nnsearch Object
#'
#' Convert a nearest neighbor search result to a sparse adjacency matrix.
#'
#' @param x An object of class "nnsearch".
#' @param idim The number of rows in the resulting matrix.
#' @param jdim The number of columns in the resulting matrix.
#' @param return_triplet Logical; whether to return triplet format.
#'
#' @return A sparse Matrix representing the adjacency matrix.
#'
#' @method adjacency nnsearch
#' @export
adjacency.nnsearch <- function(x, idim=nrow(x$idx), jdim=max(x$idx), return_triplet=FALSE) {
  indices_to_sparse(as.matrix(x$idx), as.matrix(x$dist), return_triplet=FALSE, idim=idim, jdim=jdim)
}


#' Find Nearest Neighbors Using nnsearcher
#'
#' Search for the k nearest neighbors using a pre-built nnsearcher object.
#'
#' @param x An object of class "nnsearcher".
#' @param query A matrix of query points. If NULL, searches within the original data.
#' @param k The number of nearest neighbors to find.
#' @param ... Additional arguments (currently unused).
#'
#' @return An object of class "nn_search" containing indices, distances, and labels.
#'
#' @examples
#' \donttest{
#' # Create example data and searcher
#' X <- matrix(rnorm(100), nrow=10, ncol=10)
#' searcher <- nnsearcher(X)
#' result <- find_nn(searcher, k=3)
#' }
#'
#' @method find_nn nnsearcher
#' @export
find_nn.nnsearcher <- function(x, query=NULL, k=5, ...) {
  ret <- if (!is.null(query)) {
    chk_matrix(query)
    hnsw_search(query, x$ann, k = k)
  } else {
    hnsw_search(x$X, x$ann, k = k)
  }

  ret$labels <- t(apply(ret$idx, 1,  function(i) x$labels[i]))
  search_result(x,ret)
}

#' Find Nearest Neighbors Among Subset Using nnsearcher
#'
#' Search for the k nearest neighbors within a specified subset of points.
#'
#' @param x An object of class "nnsearcher".
#' @param k The number of nearest neighbors to find.
#' @param idx A numeric vector specifying the subset of point indices to search among.
#' @param ... Additional arguments (currently unused).
#'
#' @return An object of class "nn_search" containing indices, distances, and labels.
#'
#' @method find_nn_among nnsearcher
#' @export
find_nn_among.nnsearcher <- function(x, k=5, idx, ...) {
  chk_numeric(idx)

  X1 <- x$X[idx,]
  ann <- hnsw_build(X1, x$distance, M=x$M, ef=x$ef)
  nnres <- hnsw_search(X1, ann, k=k)
  search_result(x, nnres)
}

#' Find Nearest Neighbors Among Classes
#'
#' Find the nearest neighbors within each class for a class_graph object.
#'
#' @param x A class_graph object.
#' @param X The data matrix corresponding to the graph nodes.
#' @param k The number of nearest neighbors to find.
#' @param ... Additional arguments (currently unused).
#'
#' @return A search result object containing indices, distances, and labels.
#'
#' @method find_nn_among class_graph
#' @export
find_nn_among.class_graph <- function(x, X, k=5, ...) {
  searcher <- nnsearcher(X,x$labels)
  ret <- lapply(x$class_indices, function(ind) {
    nnr <- find_nn_among.nnsearcher(searcher, k=k, ind)
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

#' Find Nearest Neighbors Between Two Sets Using nnsearcher
#'
#' Search for the k nearest neighbors from one set of points to another set.
#'
#' @param x An object of class "nnsearcher".
#' @param k The number of nearest neighbors to find.
#' @param idx1 A numeric vector specifying indices of the first set of points.
#' @param idx2 A numeric vector specifying indices of the second set of points.
#' @param restricted Logical; if TRUE, use restricted search mode.
#' @param ... Additional arguments (currently unused).
#'
#' @return An object of class "nn_search" containing indices, distances, and labels.
#'
#' @method find_nn_between nnsearcher
#' @export
find_nn_between.nnsearcher <- function(x, k=5, idx1, idx2, restricted=FALSE, ...) {
  chk_numeric(idx1)
  chk_numeric(idx2)

  ret <- if (!restricted) {
    X1 <- x$X[idx1,]
    X2 <- x$X[idx2,]
    ann <- hnsw_build(X1, x$distance, M=x$M, ef=x$ef)
    nnres <- hnsw_search(X2, ann, k=k)
    nnres$idx <- do.call(rbind, lapply(1:nrow(nnres$idx), function(i) {
      idx1[nnres$idx[i,]]
    }))
    nnres$labels <- t(apply(nnres$idx, 1,  function(i) x$labels[i]))
    nnres
  } else {
    find_nn.nnsearcher(x, x$X[idx2,])
  }

  search_result(x, ret)
}

#' Create Neighbor Graph from nnsearcher Object
#'
#' Construct a neighbor graph from nearest neighbor search results.
#'
#' @param x An object of class "nnsearcher".
#' @param query A matrix of query points. If NULL, uses original data.
#' @param k The number of nearest neighbors to find.
#' @param type The type of graph construction method.
#' @param transform The transformation method for converting distances to weights.
#' @param sigma The bandwidth parameter for the transformation.
#' @param ... Additional arguments (currently unused).
#'
#' @return A neighbor_graph object representing the constructed graph.
#'
#' @method neighbor_graph nnsearcher
#' @export
neighbor_graph.nnsearcher <- function(x, query=NULL, k=5, type=c("normal", "asym", "mutual"),
                                                           transform=c("heat", "binary", "euclidean",
                                                              "normalized", "cosine", "correlation"), sigma=1, ...) {

  nn <- find_nn.nnsearcher(x,k=(k+1))
  D <- nn$dist[,2:(k+1)]

  nnd <- sqrt(nn$dist[, 2:ncol(nn$dist),drop=FALSE])
  nni <- nn$idx[, 2:ncol(nn$idx), drop=FALSE]


  hfun <- get_neighbor_fun(transform, len=ncol(x$X), sigma=sigma)
  hval <- hfun(D)

  W <- if (k == 1) {
    indices_to_sparse(as.matrix(nni), as.matrix(hval), idim=nrow(x$X), jdim=nrow(x$X))
  } else {
    indices_to_sparse(nni, hval, idim=nrow(x$X), jdim=nrow(x$X))
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
