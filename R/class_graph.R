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

#' @export
nclasses.class_graph <- function(x) {
  length(x$levels)
}


#' @export
class_means.class_graph <- function(x, X) {
  ret <- do.call(rbind, lapply(x$class_indices, function(i) {
    colMeans(X[i,,drop=FALSE])
  }))

  row.names(ret) <- names(x$class_indices)
  ret
}


interclass_locp_scatter.class_graph <- function(x, knn_graph, X, dens=NULL) {
  if (is.null(dens)) {
    dens <- node_density(x, X)
  }

}

intraclass_locp_scatter.class_graph <- function(x, knn_graph, X,dens=NULL) {
  gi <- igraph::intersection(x, knn_graph)

  if (is.null(dens)) {
    dens <- node_density(x, X)
  }



}

#' @export
intraclass_scatter.class_graph <- function(x, X, q=1.5, mc.cores=1) {
  dens <-  intraclass_density(x, X, q, mc.cores)
  ipairs <- intraclass_pairs(x)
  ret <- lapply(1:nrow(ipairs), function(n) {
    i <- ipairs[n,1]
    j <- ipairs[n,2]
    xi <- X[i,]
    xj <- X[j,]
    dij <- sum( (xi - xj)^2)
    edij <- exp(-dij/dens[i])
    w <- 1/2 * edij * (1 + edij)
    cbind(i, j, w)
  })

  triplet <- do.call(rbind, ret)
  W <- sparseMatrix(i=triplet[,1], j=triplet[,2], x=triplet[,3])
  W <- (W + t(W))/2
}

#' @export
interclass_scatter.class_graph <- function(x, X, q=1.5, mc.cores=1) {
  dens <-  interclass_density(x, X, q, mc.cores)
  ipairs <- interclass_pairs(x)
  ret <- lapply(1:nrow(ipairs), function(n) {
    i <- ipairs[n,1]
    j <- ipairs[n,2]
    xi <- X[i,]
    xj <- X[j,]
    dij <- sum( (xi - xj)^2)
    edij <- exp(-dij/dens[i])
    w <- 1/2 * edij * (1 - edij)
    cbind(i, j, w)
  })

  triplet <- do.call(rbind, ret)
  W <- sparseMatrix(i=triplet[,1], j=triplet[,2], x=triplet[,3])
  W <- (W + t(W))/2
}


## Discriminative globality and locality preserving graph embedding for
## dimensionality reduction
## equations 11 (12)
#' @export
intraclass_density.class_graph <- function(x, X, q=1.5, mc.cores=1) {
  assert_that(nrow(X) == nvertices(x))
  labs <- x$labels


  #1/(n_c[i]^q) * sum((x_i - x_j)^2) for all x_i
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

#' @export
interclass_density.class_graph <- function(x, X, q=2, mc.cores=1) {
  assert_that(nrow(X) == nvertices(x))
  labs <- x$labels

  vind <- 1:nvertices(x)
  sapply(1:nrow(X), function(i) {
    #nabes <- neighbors(x, i)[[1]]
    #non <- vind[vind %in% nabes]
    non <- vind[-match(nabes,vind)]
    #non <- which(!(x$data[i,,drop=FALSE] == 1))
    S <- sum(sweep(X[non,], 2, X[i,], "-")^2)
    1/(length(non)^q) * S
  })

}

#' @export
intraclass_pairs.class_graph <- function(x) {
  ## convert to triplet-style matrix?
  do.call(rbind, lapply(1:length(x$class_indices), function(i) {
    expand.grid(x=x$class_indices[[i]], y=x$class_indices[[i]])
  }))

}

#' @export
interclass_pairs.class_graph <- function(x) {
  do.call(rbind, lapply(1:length(x$class_indices), function(i) {
    expand.grid(x=x$class_indices[[i]], y=unlist(x$class_indices[-i]))
  }))

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

#' @export
discriminating_distance <- function(X, k=length(labels)/2, sigma,labels) {
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

#' Compute similarity graph weighted by class structure
#'
#' @param X
#' @param k
#' @param sigma
#' @param cg
#'
#' @examples
#'
#' X <- matrix(rnorm(100*100), 100,100)
#' labels <- factor(rep(1:5, each=20))
#' cg <- class_graph(labels)
#' sigma <- .7
#'
#' @export
#'
## ref: Local similarity and diversity preserving discriminant projection for face and
## handwriting digits recognition
discriminating_simililarity <- function(X, k=length(labels)/2, sigma,cg, threshold=.01) {
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



