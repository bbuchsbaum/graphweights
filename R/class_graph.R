#' @export
#' @param labels the class label vector
#' @importFrom Matrix sparseVector tcrossprod Matrix t
#'
#' @examples
#' data(iris)
#' labels <- iris[,5]
#' cg <- class_graph(labels)
class_graph <- function(labels, sparse=TRUE) {
  labels <- as.factor(labels)
  out <- Reduce("+", lapply(levels(labels), function(lev) {
    kronecker(Matrix(labels==lev, sparse=sparse), t(Matrix(labels==lev, sparse=sparse)))
  }))

  ret <- list(
    data = out,
    labels=labels,
    class_indices=split(1:length(labels), labels),
    class_freq=table(labels),
    levels=unique(labels))

  class(ret) <- "class_graph"
  ret
}

nclasses.class_graph <- function(x) {
  length(x$levels)
}

class_means.class_graph <- function(x, X) {
  ret <- do.call(rbind, lapply(x$class_indices, function(i) {
    colMeans(X[i,,drop=FALSE])
  }))

  row.names(ret) <- names(x$class_indices)
  ret
}

intraclass_density.class_graph <- function(x, X, q=2, mc.cores=1) {
  assert_that(nrow(X) == nrow(x$data))
  labs <- x$labels
  ret <- parallel::mclapply(x$levels, function(l) {
    idx <- which(labs == l)
    xc <- X[idx,]
    sapply(1:length(idx), function(i) {
      S <- sum(sweep(xc[-i,], 2, xc[i,], "-")^2)
      1/(length(idx)^q) * S
    })
  })

  dens <- numeric(nrow(X))
  for (i in 1:nclasses(x)) {
    dens[x$class_indices[[i]]] <- ret[[i]]
  }

  dens
}


interclass_density.class_graph <- function(x, X, q=2, mc.cores=1) {
  assert_that(nrow(X) == nrow(x$data))
  labs <- x$labels

  sapply(1:nrow(X), function(i) {
    non <- which(!(x$data[i,,drop=FALSE] == 1))
    S <- sum(sweep(X[non,], 2, X[i,], "-")^2)
    1/(length(non)^q) * S
  })

}



intra_class_pairs <- function(x) {

}

inter_class_pairs <- function(x) {
}


#' @inheritParams graph_weights
#' @export
within_class_neighbors.class_graph <- function(x, X, k=1,
                                 weight_mode=c("heat", "normalized", "binary", "euclidean"),
                                 sigma,
                                 type=c("normal", "mutual", "asym"),
                                 ...) {


  type <- match.arg(type)
  weight_mode <- match.arg(weight_mode)
  assert_that(k > 0)

  if (missing(sigma)) {
    sigma <- estimate_sigma(X)
  }

  wfun <- get_neighbor_fun(weight_mode, sigma, ncol(X))

  M <- do.call(rbind, lapply(x$class_indices, function(idx) {
      k <- min(k, length(idx)-1)
      M <- weighted_knn(X[idx,], k=k, FUN=wfun, type=type,...)
      Mind <- which(M > 0, arr.ind=TRUE)
      cbind(idx[Mind[,1]], idx[Mind[,2]], M[Mind])
  }))

  W <- Matrix::sparseMatrix(i=M[,1], j=M[,2], x=M[,3], dims=c(nrow(X), nrow(X)))


}

#' @inheritParams graph_weights
#' @export
between_class_neighbors.class_graph <- function(x, X, k=1,
                                  weight_mode=c("heat", "normalized", "binary", "euclidean"),
                                  sigma,
                                  type=c("normal", "mutual", "asym"),...) {

  type <- match.arg(type)
  weight_mode <- match.arg(weight_mode)
  assert_that(k > 0)

  if (missing(sigma)) {
    sigma <- estimate_sigma(X)
  }

  wfun <- get_neighbor_fun(weight_mode, sigma, ncol(X))

  M <- do.call(rbind, lapply(x$levels, function(lev) {
    idx1 <- which(labels != lev)
    idx2 <- which(labels == lev)

    M <- weighted_knnx(X[idx1,,drop=FALSE], X[idx2,,drop=FALSE], k=k, FUN=wfun,
                       type="asym",...)
    Mind <- which(M > 0, arr.ind=TRUE)
    cbind(idx2[Mind[,1]], idx1[Mind[,2]], M[Mind])
  }))


  W <- Matrix::sparseMatrix(i=M[,1], j=M[,2], x=M[,3], dims=c(nrow(X), nrow(X)), use.last.ij = TRUE)

  if (type == "normal") {
    psparse(W, pmax)
  } else if (type == "mutual") {
    #psparse(W, pmin, return_triplet=return_triplet)
    psparse(W, pmin)
  } else if (type == "asym") {
    W
  }

}

