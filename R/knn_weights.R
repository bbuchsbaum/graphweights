
as_triplet <- function(M) {
  tm <- as(M, "dgTMatrix")
  cbind(i=tm@i+1, j=tm@j+1, x=tm@x)
}

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

#' heat_kernel
#'
#' @param x the distances
#' @param sigma the bandwidth
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


correlation_kernel <- function(x, len) {
  1 - (x^2)/(2*(len-1))
}

cosine_kernel <- function(x, sigma=1) {
  1 - (x^2)/2
}


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


#' factor_sim
#'
#' compute similarity matrix from a set of \code{factor}s in a \code{data.frame}
#'
#' @param des
#' @param method
#' @export
factor_sim <- function(des, method=c("Jaccard", "Rogers", "simple matching", "Dice")) {
  Fmat <- do.call(cbind, lapply(1:ncol(des), function(i) {
    nam <- colnames(des)[i]
    model.matrix(as.formula(paste("~ ", nam, " -1")), data=des)
  }))

  proxy::simil(Fmat, method=method)
}


weighted_factor_sim <- function(des, wts=rep(1,ncol(des))/ncol(des)) {
  wts <- wts/sum(wts)
  Fmat <- lapply(1:ncol(des), function(i) {
    nam <- colnames(des)[i]
    labs <- des[[nam]]
    label_matrix(labs, labs) * wts[i]
  })

  Reduce("+", Fmat)
}



#' @export
estimate_sigma <- function(X, prop=.25, nsamples=500) {
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

#' Convert a data matrix with n instances and p features to an n-by-n adjacency graph
#'
#' @param X the data matrix, where each row is an instance and each column is a variable. Similarity is computed over instances.
#' @param k number of neighbors (ignored when neighbor_mode is not 'epsilon')
#' @param neighbor_mode the method for assigning weights to neighbors, either "supervised", "knn", "knearest_misses", or "epsilon"
#' @param weight_mode binary (1 if neighbor, 0 otherwise),'heat', 'normalized', 'euclidean', or 'cosine'
#' @param type the nearest neighbor policy, one of: normal, mutual, asym.
#' @param sigma parameter for heat kernel \code{exp(-dist/(2*sigma^2))}
#' @param eps the neighborhood radius when neighbor_mode is `epsilon` (not implemented)
#' @param labels the class of the categories when \code{weight_mode} is \code{supervised},
#' supplied as a \code{factor} with \code{nrow(labels) == nrow(X)}
#' @export
#' @examples
#'
#' X <- matrix(rnorm(100*100), 100, 100)
#' sm <- graph_weights(X, neighbor_mode="knn",weight_mode="normalized", k=3)
#'
#' labels <- factor(rep(letters[1:4],5))
#' sm3 <- graph_weights(X, neighbor_mode="knn",k=3, labels=labels, weight_mode="cosine")
graph_weights <- function(X, k=5, neighbor_mode=c("knn", "epsilon"),
                                 weight_mode=c("heat", "normalized", "binary", "euclidean",
                                               "cosine", "correlation"),
                                 type=c("normal", "mutual", "asym"),
                                 sigma,eps=NULL, labels=NULL, ...) {

  neighbor_mode = match.arg(neighbor_mode)
  weight_mode = match.arg(weight_mode)
  type <- match.arg(type)

  if (weight_mode == "normalized" || weight_mode == "correlation") {
    X <- t(scale(t(X), center=TRUE, scale=TRUE))
  } else if (weight_mode == "cosine") {
    X <- t(apply(X, 1, function(x) x/sqrt(sum(x^2))))
  }

  if ((missing(sigma) || is.null(sigma)) && (weight_mode %in% c("heat", "normalized"))) {
    if (weight_mode == "heat") {
      sigma <- estimate_sigma(X)
    } else {
      sigma <- estimate_sigma(X)
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



#' threshold_adjacency
#'
#' extract the k-nearest neighbors from an existing adjacency matrix
#'
#' @param A adjacency matrix
#' @param k number of neighbors
#' @param type 'normal' or 'mutual' nearest neighbors
#' @param ncores number of cores to use
#' @importFrom parallel mclapply
#' @importFrom Matrix which sparseMatrix
#' @importFrom assertthat assert_that
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


# weighted_knn
#
# @param X the data matrix, where rows are instances and columns are features
# @param eps the e-radius used to threshold neighbor inclusion
# @param FUN the function used to convert euclidean distances to similarities
# @importFrom assertthat assert_that
# @importFrom FNN get.knn
# @importFrom Matrix t
# @export
# weighted_radius <- function(X, eps, FUN=heat_kernel, max_neighbour=.5*nrow(X), return_triplet=FALSE) {
#
#   if (missing(eps)) {
#     k <- round(.1 * nrow(X))
#     nn <- FNN::get.knn(X, k=k)
#     eps <- 2* min(nn$nn.dist[,2]^2)
#   }
#
#
#
#   browser()
#   res <- rflann::RadiusSearch(X, X, radius=eps, max_neighbour=max_neighbour)
#   hval <- lapply(res$distances, function(x) FUN(sqrt(x))
#
#   W <- indices_to_sparse(nn$nn.index, hval)
#
#   if (type == "normal") {
#     psparse(W, pmax, return_triplet=return_triplet)
#   } else {
#     psparse(W, pmin, return_triplet=return_triplet)
#   }
# }

weighted_knnx <- function(X, query, k=5, FUN=heat_kernel, type=c("normal", "mutual", "asym"),
                          as=c("igraph", "sparse", "index_sim")) {

  assert_that(k > 0 && k <= nrow(X))

  as <- match.arg(as)
  type <- match.arg(type)
  #nn <- FNN::get.knnx(X, query, k=k)
  #nnd <- nn$nn.dist
  #hval <- FUN(nnd)

  nn <- rflann::Neighbour(as.matrix(query), X, k=k+1)
  nnd <- sqrt(nn$distances[, 2:ncol(nn$distances),drop=FALSE] + 1e-16)
  nni <- nn$indices[, 2:ncol(nn$indices),drop=FALSE]
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

#' weighted_knn
#'
#' @param X the data matrix, where rows are instances and columns are features
#' @param k the number of nearest neighbors
#' @param FUN the function used to convert euclidean distances to similarities
#' @param type whether to compute 'normal' or 'mutual' k-nearest neighbors
#' @param as return as type \code{sparseMatrix} ("sparse") or \code{igraph} ("graph")
#' @importFrom assertthat assert_that
#' @importFrom FNN get.knn
#' @importFrom Matrix t
#' @export
#'
#' @examples
#'
#' X <- matrix(rnorm(10*10), 10,10)
#' w <- weighted_knn(X,k=5)
#'
weighted_knn <- function(X, k=5, FUN=heat_kernel,
                         type=c("normal", "mutual", "asym"),
                         as=c("igraph", "sparse"),
                         ...) {
  assert_that(k > 0 && k <= nrow(X))
  as <- match.arg(as)


  type <- match.arg(type)
  #nn <- FNN::get.knn(X, k=k)
  #nn <- nabor::knn(X, k=k)
  nn <- rflann::Neighbour(X, X,k=k+1, ...)
  #nnd <- nn$nn.dist + 1e-16

  nnd <- sqrt(nn$distances[, 2:ncol(nn$distances),drop=FALSE] + 1e-16)
  nni <- nn$indices[, 2:ncol(nn$indices), drop=FALSE]
  hval <- FUN(nnd)

  #W <- indices_to_sparse(nn$nn.index, hval)
  W <- if (k == 1) {
    indices_to_sparse(as.matrix(nni), as.matrix(hval), idim=nrow(X), jdim=nrow(X))
  } else {
    indices_to_sparse(nni, hval, idim=nrow(X), jdim=nrow(X))
  }

  ##bW <- W != 0
  ##tmp <- as(x$G, "dgTMatrix")
  ## should return the raw distances?

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

#' psparse
#'
#' apply a function (e.g. max) to each pair of nonzero elements in a \code{sparseMatrix}
#'
#' @param M
#' @param FUN
#' @param return_triplet
#' @importFrom Matrix which
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



# psparse <- function(..., fun=c("max", "min"), na.rm=FALSE, return_triplet=FALSE) {
#   fun <- match.arg(fun)
#   # check that all matrices have conforming sizes
#   num.rows <- unique(sapply(list(...), nrow))
#   num.cols <- unique(sapply(list(...), ncol))
#   stopifnot(length(num.rows) == 1)
#   stopifnot(length(num.cols) == 1)
#
#   cat.summary <- do.call(rbind, lapply(list(...), summary))
#
#
#   out.summary <- if (fun == "min") {
#     aggregate(x ~ i + j, data = cat.summary, FUN=function(x) {
#       if (length(x) == 1) 0 else x[1]
#     })
#   } else {
#     aggregate(x ~ i + j, data = cat.summary, FUN=max, na.rm)
#   }
#
#   if (return_triplet) {
#     cbind(i=out.summary$i, j=out.summary$j, x=out.summary$x)
#   } else {
#     sparseMatrix(i = out.summary$i,
#                  j = out.summary$j,
#                  x = out.summary$x,
#                  dims = c(num.rows, num.cols))
#   }
# }
#
# pmin.sparse <- function(..., na.rm = FALSE, return_triplet=FALSE) { psparse(..., fun="min", return_triplet=return_triplet) }
#
# pmax.sparse <- function(..., na.rm = FALSE, return_triplet=FALSE) { psparse(..., fun="max", return_triplet=return_triplet ) }
#

