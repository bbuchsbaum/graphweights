#' Compute Repulsion Weight Between Feature Vectors
#'
#' @description
#' Computes a repulsion weight between two feature vectors, where the weight
#' decreases as the distance between vectors increases, normalized by their magnitudes.
#'
#' @details
#' The repulsion weight is computed as:
#' w = 1 / (sigma + ||x1 - x2||^2 / (||x1|| + ||x2||))
#' where:
#' * ||x1 - x2|| is the Euclidean distance between vectors
#' * ||x1||, ||x2|| are the L2 norms of the vectors
#' * sigma is a smoothing parameter
#'
#' This weight function has the following properties:
#' * Decreases with increasing distance between vectors
#' * Normalized by vector magnitudes to reduce scale sensitivity
#' * Bounded between 0 and 1/(sigma)
#'
#' @param x1 Numeric vector representing the first feature vector
#' @param x2 Numeric vector representing the second feature vector
#' @param sigma Numeric smoothing parameter (default: 10)
#'
#' @return A numeric value between 0 and 1/sigma representing the repulsion weight
#'
#' @keywords internal
repulse_weight <- function(x1, x2, sigma=10) {
  num <- sum((x1 - x2)^2)
  denom <- sqrt(sum(x1^2)) + sqrt(sum(x2^2))
  1/(sigma + num/denom)
}

#' Create a Repulsion Graph
#'
#' @description
#' Creates a repulsion graph that emphasizes edges between vertices of different
#' classes while suppressing edges within the same class. This is particularly
#' useful for semi-supervised learning and graph-based classification tasks where
#' class boundaries need to be emphasized.
#'
#' @details
#' A repulsion graph modifies an input similarity graph by incorporating class
#' information to enhance between-class edges and suppress within-class edges.
#' This is achieved by:
#'
#' 1. Binary method:
#'    * Removes all edges between vertices of the same class
#'    * Preserves edges between different classes
#'
#' 2. Weighted method:
#'    * Scales edge weights by a normalization factor
#'    * Removes within-class edges
#'    * Preserves between-class edge structure
#'
#' Common use cases:
#' * Semi-supervised learning with graph-based methods
#' * Enhancing class separation in spectral clustering
#' * Feature selection based on class discrimination
#' * Graph-based dimensionality reduction
#'
#' @param W A weighted adjacency matrix (class Matrix) or neighbor_graph object
#'   representing the base similarity graph
#' @param cg A class_graph object indicating pairwise relationships among node classes
#' @param method Character string specifying the repulsion method:
#'   * "binary": Creates a binary repulsion graph
#'   * "weighted": Creates a weighted repulsion graph (default)
#' @param threshold Numeric threshold for filtering edge weights (default: 0).
#'   Edges with weights <= threshold are removed
#' @param norm_fac Numeric normalization factor for edge weights in weighted method
#'   (default: 1). Higher values reduce edge weights
#'
#' @return A neighbor_graph object representing the repulsion graph where:
#'   * Edges between same-class vertices are removed
#'   * Between-class edges are preserved (binary) or scaled (weighted)
#'   * The graph maintains the same number of vertices as the input
#'
#' @examples
#' # Example 1: Using iris dataset
#' data(iris)
#' X <- as.matrix(iris[, 1:4])
#' labels <- iris$Species
#'
#' # Create base similarity graph
#' W <- graph_weights(X, k = 6)
#'
#' # Create class graph
#' cg <- class_graph(labels)
#'
#' # Create binary repulsion graph
#' R_bin <- repulsion_graph(W, cg, method = "binary")
#'
#' # Create weighted repulsion graph
#' R_wt <- repulsion_graph(W, cg, method = "weighted", norm_fac = 2)
#'
#' # Example 2: Custom data with multiple classes
#' set.seed(123)
#' X <- matrix(rnorm(100 * 10), 100, 10)  # 100 samples, 10 features
#' labels <- factor(rep(1:5, each = 20))   # 5 classes
#'
#' # Create graphs
#' W <- graph_weights(X, k = 8)
#' cg <- class_graph(labels)
#' R <- repulsion_graph(W, cg, method = "weighted", threshold = 0.1)
#'
#' @seealso
#' * \code{\link{class_graph}} for creating class relationship graphs
#' * \code{\link{graph_weights}} for creating base similarity graphs
#' * \code{\link{neighbor_graph}} for the base graph class
#'
#' @importFrom assertthat assert_that
#' @export
repulsion_graph <- function(W, cg, method=c("binary", "weighted"), threshold=0,
                          norm_fac=1) {
  method <- match.arg(method)

  # Input validation and conversion
  if (inherits(W, "neighbor_graph")) {
    W <- adjacency(W)
    W <- W * (W > threshold)
  } else if (inherits(W, "Matrix")) {
    W <- W * (W > threshold)
  } else {
    stop("`W` must be a `Matrix` or `neighbor_graph` object")
  }

  # Create repulsion graph based on method
  res <- if (method == "binary") {
    neighbor_graph((W * !adjacency(cg)) > 0)
  } else {
    R <- W * !adjacency(cg)
    ind <- which(R != 0, arr.ind=TRUE)
    R[ind] <- W[ind]/norm_fac
    neighbor_graph(R)
  }

  res
}
