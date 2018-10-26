

#' internal function
indices_to_sparse <- function(nn.index, hval, return_triplet=FALSE, idim=nrow(nn.index), jdim=nrow(nn.index)) {
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

cosine_kernel <- function(x, sigma=1) {
  1 - (x^2)/2
}


#' factor_sim
#'
#' compute similarity matrix from a set of \code{factor}s in a \code{data.frame}
#'
#' @param des
#' @param method
#' @export
factor_sim <- function(des, method=c("Jaccard", "Rogers", "simple matching", "Dice")) {
  Fmat <- do.call(cbind, lapply(names(des), function(nam) {
    model.matrix(as.formula(paste("~ ", nam, " -1")), data=des)
  }))

  proxy::simil(Fmat, method=method)
}


#'
#' @export
discriminating_simililarity <- function(X, k=length(labels)/2, sigma,labels) {
  Wknn <- graph_weights(X)

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
  hw <- heat_kernel(Wall[wind], sigma)
  hb <- heat_kernel(Wall[bind], sigma)

  Wall[wind] <- hw * (1+hw)
  Wall[bind] <- hb * (1-hb)
  Wall
}


#' @inheritParams graph_weights
#' @export
within_class_weights <- function(X, k=1, labels, weight_mode=c("heat", "normalized", "binary", "euclidean"),...) {
  labels <- as.factor(labels)
  Sw <- graph_weights(X, k=k, weight_mode=weight_mode, neighbor_mode="supervised", labels=labels,...)
}

#' @export
between_class_weights <- function(X, k=1, labels, weight_mode=c("heat", "normalized", "binary", "euclidean"), ...) {
  labels <- as.factor(labels)
  Sb <- graph_weights(X, k=k, weight_mode=weight_mode, neighbor_mode="knearest_misses", labels=labels,...)
}

#' @export
estimate_sigma <- function(X, prop=.1, nsamples=500) {
  ## estimate sigma
  if (nrow(X) <= 500) {
    nsamples <- nrow(X)
    sam <- 1:nrow(X)
  } else {
    nsamples <- min(nsamples, nrow(X))
    sam <- sample(1:nrow(X), nsamples)
  }


  d <- dist(X[sam,])
  quantile(d[d!=0],prop)

}

#' Convert a data matrix with n instances and p features to an n-by-n adjacency matrix
#'
#' @param X the data matrix, where each row is an instance and each column is a variable. Similarity is computed over instances.
#' @param k number of neighbors
#' @param neighbor_mode the method for assigning weights to neighbors, either "supervised", "knn", "knearest_misses", or "epsilon"
#' @param weight_mode binary (1 if neighbor, 0 otherwise),'heat', 'normalized', 'euclidean', or 'cosine'
#' @param sigma parameter for heat kernel \code{exp(-dist/(2*sigma^2))}
#' @param labels the class of the categories when \code{weight_mode} is \code{supervised}, supplied as a \code{factor} with \code{nrow(labels) == nrow(X)}
#' @export
#' @examples
#'
#' X <- matrix(rnorm(20*10), 20, 10)
#' sm <- graph_weights(X, neighbor_mode="knn",k=3)
#'
#' labels <- factor(rep(letters[1:4],5))
#' sm2 <- graph_weights(X, neighbor_mode="supervised",k=3, labels=labels)
#'
#' sm3 <- graph_weights(X, neighbor_mode="knearest_misses",k=3, labels=labels, weight_mode="binary")
graph_weights <- function(X, k=5, neighbor_mode=c("knn", "supervised", "knearest_misses", "epsilon"),
                                 weight_mode=c("heat", "normalized", "binary", "euclidean", "cosine"),
                                 type=c("normal", "mutual", "asym"),
                                 sigma,eps=NULL, labels=NULL) {

  neighbor_mode = match.arg(neighbor_mode)
  weight_mode = match.arg(weight_mode)
  type <- match.arg(type)


  if (weight_mode == "normalized") {
    X <- t(scale(t(X), center=TRUE, scale=TRUE))
  } else if (weight_mode == "cosine") {
    X <- t(apply(X, 1, function(x) x/sqrt(sum(x^2))))
  }

  if (missing(sigma) || is.null(sigma)) {
    sigma <- estimate_sigma(X)
    message("sigma is ", sigma)
  }

  wfun <- if (weight_mode == "heat") {
    function(x) heat_kernel(x, sigma)
  } else if (weight_mode == "binary") {
    function(x) (x>0) * 1
  } else if (weight_mode == "normalized") {
    function(x) normalized_heat_kernel(x,sigma,ncol(X))
  } else if (weight_mode == "euclidean") {
    function(x) x
  } else if (weight_mode == "cosine") {
    function(x) cosine_kernel(x)
  }

  if (neighbor_mode == "knn") {
    W <- weighted_knn(X, k, FUN=wfun, type=type)
  } else if (neighbor_mode == "supervised") {
    assertthat::assert_that(!is.null(labels))
    labels <- as.factor(labels)
    if (k == 0) {
      W <- label_matrix2(labels, labels)
    } else if (k > 0) {

      M <- do.call(rbind, lapply(levels(labels), function(lev) {
        idx <- which(labels == lev)
        k <- min(k, length(idx)-1)
        M <- weighted_knn(X[idx,], k=k, FUN=wfun, type=type)
        Mind <- which(M > 0, arr.ind=TRUE)
        cbind(idx[Mind[,1]], idx[Mind[,2]], M[Mind])
      }))

      Matrix::sparseMatrix(i=M[,1], j=M[,2], x=M[,3], dims=c(nrow(X), nrow(X)))
    }

  } else if (neighbor_mode == "knearest_misses") {
    assertthat::assert_that(!is.null(labels) && k > 0,
                            msg="when 'neighbor_mode' is 'knearest_misses',
                            k must be greater than 0 and labels cannot be NULL")
    labels <- as.factor(labels)
    M <- do.call(rbind, lapply(levels(labels), function(lev) {
      idx1 <- which(labels != lev)
      idx2 <- which(labels == lev)

      M <- weighted_knnx(X[idx1,,drop=FALSE], X[idx2,,drop=FALSE], k=k, FUN=wfun, type="asym")
      Mind <- which(M > 0, arr.ind=TRUE)
      cbind(idx2[Mind[,1]], idx1[Mind[,2]], M[Mind])
    }))


    W <- Matrix::sparseMatrix(i=M[,1], j=M[,2], x=M[,3], dims=c(nrow(X), nrow(X)), use.last.ij = TRUE)

    if (type == "normal") {
      psparse(W, pmax)
    } else if (type == "mutual") {

      #psparse(W, pmax, return_triplet=return_triplet)
      psparse(W, pmax)
    } else if (type == "mutual") {
      #psparse(W, pmin, return_triplet=return_triplet)
      psparse(W, pmin)
    } else if (type == "asym") {
      W
    }

  } else if (neighbor_mode == "epsilon") {
    stop("epsilon not implemented")
  }
}



#' sim_from_adj
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
sim_from_adj <- function(A, k=5, type=c("normal", "mutual"), ncores=1) {
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

weighted_knnx <- function(X, query, k=5, FUN=heat_kernel, type=c("normal", "mutual", "asym"), return_triplet=FALSE) {
  assert_that(k > 0 && k <= nrow(X))

  type <- match.arg(type)
  nn <- FNN::get.knnx(X, query, k=k)

  nnd <- nn$nn.dist

  hval <- FUN(nnd)

  W <- indices_to_sparse(nn$nn.index, hval, idim=nrow(nn$nn.index), jdim=nrow(X), return_triplet=return_triplet)

  if (type == "normal") {
    psparse(W, pmax, return_triplet=return_triplet)
  } else if (type == "mutual") {
    psparse(W, pmin, return_triplet=return_triplet)
  } else if (type == "asym") {
    W
  }
}

#' weighted_knn
#'
#' @param X the data matrix, where rows are instances and columns are features
#' @param k the number of nearest neighbors
#' @param FUN the function used to convert euclidean distances to similarities
#' @param type whether to compute 'normal' or 'mutual' k-nearest neighbors
#' @importFrom assertthat assert_that
#' @importFrom FNN get.knn
#' @importFrom Matrix t
#' @export
weighted_knn <- function(X, k=5, FUN=heat_kernel, type=c("normal", "mutual", "asym"), return_triplet=FALSE) {
  assert_that(k > 0 && k <= nrow(X))

  type <- match.arg(type)
  nn <- FNN::get.knn(X, k=k)

  nnd <- nn$nn.dist + 1e-16

  hval <- FUN(nnd)

  W <- indices_to_sparse(nn$nn.index, hval)

  if (type == "normal") {
    psparse(W, pmax, return_triplet=return_triplet)
  } else if (type == "mutual") {
    psparse(W, pmin, return_triplet=return_triplet)
  } else if (type == "asym") {
    W
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
    sm <- sparseMatrix(i=c(ind[,1],ind[,2]), j=c(ind[,2],ind[,1]), x=rep(FUN(x1,x2),2), dims=dim(M), use.last.ij=TRUE)
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

