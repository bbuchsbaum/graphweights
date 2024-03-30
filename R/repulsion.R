
#' @keywords internal
repulse_weight <- function(x1,x2, sigma=10) {
  #num <- sum((x1 - x2)^2)
  num <- sum((x1 - x2)^2)
  denom <- sqrt(sum(x1^2)) + sqrt(sum(x2^2))

  1/(sigma + num/denom)
}


#' This function creates a repulsion graph, which is a graph where edges between vertices of different classes are emphasized, while edges within the same class are suppressed.
#'
#' @param W The adjacency matrix (weighted or binary) representing the relationships between vertices in the input graph.
#' @param cg The class graph indicating the pairwise relationships among node classes.
#' @param method A character string indicating the method used to create the repulsion graph. Can be either "binary" or "weighted". Default is "binary".
#' @param threshold A numeric value used to filter the adjacency matrix. Edges with weights below this threshold will be ignored. Default is 0.
#' @param norm_fac A numeric value used to normalize the weights in the adjacency matrix when the method is set to "weighted". Default is 1.
#'
#' @return A `neighbor_graph` object representing the repulsion graph.
#'
#' @examples
#' # Example with random data
#' X <- matrix(rnorm(100*100), 100, 100)
#' labels <- factor(rep(1:5, each=20))
#' cg <- class_graph(labels)
#' W <- graph_weights(X, k=10)
#' R <- repulsion_graph(W, cg, method="weighted")
#'
#' # Example with iris dataset
#' data(iris)
#' X2 <- as.matrix(iris[, 1:4])
#' labels <- iris[, 5]
#' cg <- class_graph(labels)
#' W <- graph_weights(X2, k=6)
#' R <- repulsion_graph(W, cg, method="weighted")
#'
#' @export
repulsion_graph <- function(W, cg, method=c("binary", "weighted"), threshold=0,
                            norm_fac=1) {
  method <- match.arg(method)
  #cg <- class_graph(labels)
  #acg <- 1- (adjacency(cg) > threshold)

  if (inherits(W, "neighbor_graph")) {
    W <- adjacency(W)
    W <- W * (W > threshold)
  } else if (inherits(W, "Matrix")) {
    W <- W * (W > threshold)
  } else {
    stop("`W` must be a `Matrix` or `neighbor_graph` object")
  }

  res <- if (method == "binary") {
    neighbor_graph( (W * !adjacency(cg)) > 0)
    #W * !adjacency(cg)
  } else {
    #assert_that(!is.null(X))
    #assert_that(nrow(X) == nrow(W))
    R <- W * !adjacency(cg)
    ind <- which(R != 0, arr.ind=TRUE)
    #R <- acg * (W != 0)
    R[ind] <- W[ind]/norm_fac
    neighbor_graph(R)
  }

}
