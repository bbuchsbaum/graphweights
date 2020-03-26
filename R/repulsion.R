
repulse_weight <- function(x1,x2, sigma=10) {
  #num <- sum((x1 - x2)^2)
  num <- sum((x1 - x2)^2)
  denom <- sqrt(sum(x1^2)) + sqrt(sum(x2^2))

  1/(sigma + num/denom)
}

#' construct a "repulsion graph"
#'
#' @param labels a set of labels defining the classes
#' @param W the neighborhood graph (weighted or binary)
#' @param X the original data matrix where \code{nrow(X) == nrow(W)}
#' @param sigma bandwidth of the similarity function
#' @export
repulsion_graph <- function(labels, W, X=NULL, sigma=10, method=c("binary", "weighted")) {
  method <- match.arg(method)
  cg <- class_graph(labels)
  acg <- 1-cg

  res <- if (type == "binary") {
    acg * (W != 0)
  } else {
    assert_that(!is.null(X))
    assert_that(nrow(X) == nrow(W))
    ind <- which((acg * W) != 0, arr.ind=TRUE)
    R <- acg * (W != 0)
    wts <- sapply(1:nrow(ind), function(i) {
      repulse_weight(X[ind[i,1],], X[ind[i,2],], sigma)
    })
  }

}
