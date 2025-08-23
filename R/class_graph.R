#' Construct a Class Graph
#'
#' A graph in which members of the same class have edges.
#'
#' @param labels A vector of class labels.
#' @param sparse A logical value, indicating whether to use sparse matrices in the computation. Default is TRUE.
#'
#' @return A class_graph object, which is a list containing the following components:
#' \itemize{
#'   \item{adjacency}{A matrix representing the adjacency of the graph.}
#'   \item{params}{A list of parameters used in the construction of the graph.}
#'   \item{labels}{A vector of class labels.}
#'   \item{class_indices}{A list of vectors, each containing the indices of elements belonging to a specific class.}
#'   \item{class_freq}{A table of frequencies for each class.}
#'   \item{levels}{A vector of unique class labels.}
#'   \item{classes}{A character string indicating the type of graph ("class_graph").}
#' }
#'
#' @importFrom Matrix sparseVector tcrossprod Matrix t
#'
#' @examples
#' data(iris)
#' labels <- iris[,5]
#' cg <- class_graph(labels)
#'
#' @export
class_graph <- function(labels, sparse=TRUE) {
  labels <- as.factor(labels)
  out <- Reduce("+", lapply(levels(labels), function(lev) {
    kronecker(Matrix(labels==lev, sparse=sparse), t(Matrix(labels==lev, sparse=sparse)))
  }))

  ret <- neighbor_graph(
    out,
    params = list(weight_mode="binary", neighbor_mode="supervised"),
    labels=labels,
    class_indices=split(1:length(labels), labels),
    class_freq=table(labels),
    levels=unique(labels),
    classes="class_graph"
  )

  ret
}



#' Number of Classes for class_graph Objects
#'
#' Compute the number of classes in a class_graph object.
#'
#' @param x A class_graph object.
#'
#' @return The number of classes in the class_graph.
#'
#' @method nclasses class_graph
#' @export
nclasses.class_graph <- function(x) {
  length(x$levels)
}


#' Class Means for class_graph Objects
#'
#' Compute the mean of each class for a class_graph object.
#'
#' @param x A class_graph object.
#' @param X The data matrix corresponding to the graph nodes.
#' @param ... Additional arguments (currently ignored).
#'
#' @return A matrix where each row represents the mean values for each class.
#'
#' @method class_means class_graph
#' @export
class_means.class_graph <- function(x, X, ...) {
  ret <- do.call(rbind, lapply(x$class_indices, function(i) {
    colMeans(X[i,,drop=FALSE])
  }))

  row.names(ret) <- names(x$class_indices)
  ret
}


#' Heterogeneous Neighbors for class_graph Objects
#'
#' Compute the neighbors between different classes for a class_graph object.
#'
#' @param x A class_graph object.
#' @param X The data matrix corresponding to the graph nodes.
#' @param k The number of nearest neighbors to find.
#' @param weight_mode Method for weighting edges (e.g., "heat", "binary", "euclidean").
#' @param sigma Scaling factor for heat kernel if `weight_mode="heat"`.
#' @param ... Additional arguments passed to weight function.
#'
#' @return A neighbor_graph object representing the between-class neighbors.
#' @importFrom FNN get.knnx
#' @importFrom Matrix sparseMatrix
#' @export
heterogeneous_neighbors <- function(x, X, k, weight_mode="binary", sigma=1, ...) {
  all_indices <- 1:nrow(X)
  class_indices <- x$class_indices
  all_triplets <- list()

  weight_fun <- get_neighbor_fun(weight_mode, sigma=sigma, ...)

  for (lev in names(class_indices)) {
    query_idx <- class_indices[[lev]]
    # Ensure data_idx contains valid indices within the bounds of X
    data_idx <- setdiff(all_indices, query_idx)
    data_idx <- data_idx[data_idx >= 1 & data_idx <= nrow(X)]
    query_idx <- query_idx[query_idx >= 1 & query_idx <= nrow(X)]

    if (length(data_idx) == 0 || length(query_idx) == 0) next

    # Ensure k is not larger than the number of available data points
    actual_k <- min(k, length(data_idx))
    if (actual_k == 0) next

    knn_result <- FNN::get.knnx(X[data_idx, , drop = FALSE],
                                 X[query_idx, , drop = FALSE],
                                 k = actual_k)

    # knn_result$nn.index gives indices relative to data_idx
    # knn_result$nn.dist gives distances

    weights <- weight_fun(knn_result$nn.dist)

    # Prepare triplets (i, j, weight)
    row_indices_orig <- rep(query_idx, actual_k)
    # Map nn.index back to original indices
    col_indices_orig <- data_idx[as.vector(knn_result$nn.index)]

    triplets <- data.frame(
      i = row_indices_orig,
      j = col_indices_orig,
      weight = as.vector(weights)
    )
    all_triplets[[lev]] <- triplets
  }

  if (length(all_triplets) == 0) {
    # Return an empty graph if no neighbors found
    adj <- Matrix::sparseMatrix(i={}, j={}, dims=c(nrow(X), nrow(X)))
  } else {
    final_triplets <- do.call(rbind, all_triplets)
    # Ensure symmetry by adding transposed triplets and averaging/maxing weights if needed
    # Or just build directed and make symmetric later
    adj <- Matrix::sparseMatrix(i = final_triplets$i,
                                j = final_triplets$j,
                                x = final_triplets$weight,
                                dims = c(nrow(X), nrow(X)))
     # Make symmetric - choose max weight for shared edges
    adj <- pmax(adj, Matrix::t(adj))
  }

  neighbor_graph(adj, params = list(weight_mode=weight_mode, neighbor_mode="heterogeneous", k=k, sigma=sigma))
}

#' Homogeneous Neighbors for class_graph Objects
#'
#' Compute the neighbors within the same class for a class_graph object.
#'
#' @param x A class_graph object.
#' @param X The data matrix corresponding to the graph nodes.
#' @param k The number of nearest neighbors to find.
#' @param weight_mode Method for weighting edges (e.g., "heat", "binary", "euclidean").
#' @param sigma Scaling factor for heat kernel if `weight_mode="heat"`.
#' @param ... Additional arguments passed to weight function.
#'
#' @return A neighbor_graph object representing the within-class neighbors.
#' @importFrom FNN get.knn
#' @importFrom Matrix sparseMatrix
#' @export
homogeneous_neighbors <- function(x, X, k, weight_mode="heat", sigma=1, ...) {
  class_indices <- x$class_indices
  all_triplets <- list()
  N <- nrow(X)

  weight_fun <- get_neighbor_fun(weight_mode, sigma=sigma, ...)

  for (lev in names(class_indices)) {
    idx_subset <- class_indices[[lev]]
    n_subset <- length(idx_subset)

    # Skip if class is too small to find k neighbors
    if (n_subset <= k) next

    knn_result <- FNN::get.knn(X[idx_subset, , drop = FALSE], k = k)

    # knn_result$nn.index gives indices relative to idx_subset
    # knn_result$nn.dist gives distances

    weights <- weight_fun(knn_result$nn.dist)

    # Row indices within the subset (repeated k times)
    row_indices_subset <- rep(1:n_subset, k)
    # Column indices within the subset (from nn.index)
    col_indices_subset <- as.vector(knn_result$nn.index)

    # Map back to original indices
    row_indices_orig <- idx_subset[row_indices_subset]
    col_indices_orig <- idx_subset[col_indices_subset]

    triplets <- data.frame(
      i = row_indices_orig,
      j = col_indices_orig,
      weight = as.vector(weights)
    )
    all_triplets[[lev]] <- triplets
  }

  if (length(all_triplets) == 0) {
    # Return an empty graph
    adj <- Matrix::sparseMatrix(i={}, j={}, dims=c(N, N))
  } else {
    final_triplets <- do.call(rbind, all_triplets)

    # Build sparse matrix from triplets
    adj <- Matrix::sparseMatrix(i = final_triplets$i,
                                j = final_triplets$j,
                                x = final_triplets$weight,
                                dims = c(N, N))
    # Make symmetric: take the max weight if edge (i,j) and (j,i) both exist
    adj <- pmax(adj, Matrix::t(adj))
  }

  neighbor_graph(adj, params = list(weight_mode=weight_mode, neighbor_mode="homogeneous", k=k, sigma=sigma))
}


#' Within-Class Neighbors for class_graph Objects
#'
#' Compute the within-class neighbors of a class_graph object.
#'
#' @param x A class_graph object.
#' @param ng A neighbor graph object.
#' @param ... Additional arguments (currently ignored).
#'
#' @return A neighbor_graph object representing the within-class neighbors of the input class_graph.
#'
#' @method within_class_neighbors class_graph
#' @export
within_class_neighbors.class_graph <- function(x, ng, ...) {
  Ac <- adjacency(x)
  An <- adjacency(ng)
  Aout <- Ac * An
  neighbor_graph(Aout)
}


#' Between-Class Neighbors for class_graph Objects
#'
#' Compute the between-class neighbors of a class_graph object.
#'
#' @param x A class_graph object.
#' @param ng A neighbor_graph object.
#' @param ... Additional arguments (currently ignored).
#'
#' @return A neighbor_graph object representing the between-class neighbors.
#'
#' @method between_class_neighbors class_graph
#' @export
between_class_neighbors.class_graph <- function(x, ng, ...) {
  Ac <- !adjacency(x)
  An <- adjacency(ng)
  Aout <- Ac * An
  neighbor_graph(Aout)
}


#' Compute Discriminating Distance for Similarity Graph
#'
#' This function computes a discriminating distance matrix for the similarity graph based on the class labels.
#' It adjusts the similarity graph by modifying the weights within and between classes, making it more suitable for
#' tasks like classification and clustering.
#'
#' @param X A numeric matrix or data frame containing the data points.
#' @param k An integer representing the number of nearest neighbors to consider.
#' @param sigma A numeric value representing the scaling factor for the heat kernel. If not provided, it will be estimated.
#' @param labels A factor or numeric vector containing the class labels for each data point.
#'
#' @return A discriminating distance matrix in the form of a numeric matrix.
#'
#' @examples
#' X <- matrix(rnorm(100*100), 100, 100)
#' labels <- factor(rep(1:5, each=20))
#' sigma <- 0.7
#' D <- discriminating_distance(X, k=length(labels)/2, sigma, labels)
#'
#' @export
discriminating_distance <- function(X, k=length(labels)/2, sigma, labels) {
  #Wknn <- graph_weights(X)

  if (missing(sigma)) {
    sigma <- estimate_sigma(X, prop=.1)
  }

  Wall <- graph_weights(X, k=k, weight_mode="euclidean", neighbor_mode="knn")

  Ww <- label_matrix2(labels, labels)
  Wb <- label_matrix2(labels, labels, type="d")

  Ww2 <- Wall * Ww
  Wb2 <- Wall * Wb

  wind <- which(Ww2 >0)
  bind <- which(Wb2 >0)

  hw <- inverse_heat_kernel(Wall[wind], sigma)
  hb <- inverse_heat_kernel(Wall[bind], sigma)

  Wall[wind] <- hw * (1-hw)
  Wall[bind] <- hb * (1+hb)
  Wall
}

#' Compute Similarity Graph Weighted by Class Structure
#'
#' This function computes a similarity graph that is weighted by the class structure of the data.
#' It is useful for preserving the local similarity and diversity within the data, making it
#' suitable for tasks like face and handwriting digits recognition.
#'
#' @param X A numeric matrix or data frame containing the data points.
#' @param k An integer representing the number of nearest neighbors to consider.
#' @param sigma A numeric value representing the scaling factor for the heat kernel.
#' @param cg A class_graph object computed from the labels.
#' @param threshold A numeric value representing the threshold for the class graph. Default is 0.01.
#'
#' @return A weighted similarity graph in the form of a matrix.
#'
#' @examples
#' X <- matrix(rnorm(100*100), 100, 100)
#' labels <- factor(rep(1:5, each=20))
#' cg <- class_graph(labels)
#' sigma <- 0.7
#' W <- discriminating_simililarity(X, k=length(labels)/2, sigma, cg)
#'
#' @references
#' Local similarity and diversity preserving discriminant projection for face and
#' handwriting digits recognition
#'
#' @export
discriminating_simililarity <- function(X, k, sigma, cg, threshold=.01) {
  #Wknn <- graph_weights(X)

  Wall <- graph_weights(X, k=k, weight_mode="heat", neighbor_mode="knn",sigma=sigma)
  Ww <- cg
  Wb <- 1- (cg > threshold)

  Ww2 <- Wall * Ww
  Wb2 <- Wall * Wb

  wind <- which(Ww2 >0)
  bind <- which(Wb2 >0)

  hw <- heat_kernel(Wall[wind], sigma)
  hb <- heat_kernel(Wall[bind], sigma)

  Wall[wind] <- hw * (1+hw)
  Wall[bind] <- hb * (1-hb)
  Wall
}



