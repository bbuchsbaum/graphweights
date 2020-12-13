
repulse_weight <- function(x1,x2, sigma=10) {
  #num <- sum((x1 - x2)^2)
  num <- sum((x1 - x2)^2)
  denom <- sqrt(sum(x1^2)) + sqrt(sum(x2^2))

  1/(sigma + num/denom)
}

#' construct a "repulsion graph"
#'
#' @param W the adjacency matrix (weighted or binary)
#' @param cg the class graph indicating the pairwise relationships among node classes
#' @param threshold
#' @export
#' @examples
#'
#' X <- matrix(rnorm(100*100), 100,100)
#' labels <- factor(rep(1:5, each=20))
#' cg <- class_graph(labels)
#' W <- graph_weights(X, k=10)
#' R <- repulsion_graph(W, cg, method="weighted")
repulsion_graph <- function(W, cg, method=c("binary", "weighted"), threshold=.001, norm_fac=10) {
  method <- match.arg(method)
  #cg <- class_graph(labels)
  acg <- 1- (cg > threshold)

  if (inherits(W, "neighbor_graph")) {
    W <- weight_matrix(W)
  }

  res <- if (method == "binary") {
    acg * (W != 0)
  } else {
    #assert_that(!is.null(X))
    #assert_that(nrow(X) == nrow(W))
    ind <- which((acg * W) != 0, arr.ind=TRUE)
    R <- acg * (W != 0)
    R[ind] <- W[ind]/norm_fac
    R
  }

}
