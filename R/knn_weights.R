

#' internal function
indices_to_sparse <- function(nn.index, hval, return_triplet=FALSE) {
  M <- do.call(rbind, lapply(1:nrow(nn.index), function(i) {
    cbind(i, nn.index[i,], hval[i,])
  }))

  if (return_triplet) {
    M
  } else {
    Matrix::sparseMatrix(i=M[,1], j=M[,2], x=M[,3], dims=c(nrow(nn.index), nrow(nn.index)))
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
normalized_heat_kernel <- function(x, sigma=1, len) {
  norm_dist <- (x^2)/(2*len)
  exp(-norm_dist/(2*sigma^2))
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


#' Convert a data matrix to a sparse similarity matrix
#'
#' @param X the data matrix, where each row is an instance and each column is a variable. Similarity is compute over instances.
#' @param neighbor_mode the method for assigning weights to neighbors, either "supervised", "knn", "knearest_misses", or "epsilon"
#' @param weight_mode binary (1 if neighbor, 0 otherwise), heat kernel, or normalized heat kernel
#' @param k number of neighbors
#' @param sigma parameter for heat kernel \code{exp(-dist/(2*sigma^2))}
#' @param labels the class of the categories when \code{weight_mode} is \code{supervised}, supplied as a \code{factor} with \code{nrow(labels) == nrow(X)}
#' @export
#' @examples
#'
#' X <- matrix(rnorm(20*10), 20, 10)
#' sm <- similarity_matrix(X, neighbor_mode="knn",k=3)
#'
#' labels <- factor(rep(letters[1:4],5))
#' sm2 <- similarity_matrix(X, neighbor_mode="supervised",k=3, labels=labels)
#'
#' sm3 <- similarity_matrix(X, neighbor_mode="knearest_misses",k=3, labels=labels, weight_mode="binary")
similarity_matrix <- function(X, neighbor_mode=c("knn", "supervised", "knearest_misses", "epsilon"),
                                 weight_mode=c("heat", "normalized", "binary"),
                                 k=5, sigma,eps=NULL, labels=NULL) {

  neighbor_mode = match.arg(neighbor_mode)
  weight_mode = match.arg(weight_mode)


  if (weight_mode == "normalized") {
    X <- t(scale(t(X), center=TRUE, scale=TRUE))
  }

  if (missing(sigma)) {
    ## estimate sigma
    if (nrow(X) <= 100) {
      nsamples <- nrow(X)
    } else {
      nsamples <- max(min(100, .25*nrow(X)), 100)

    }

    sam <- sample(1:nrow(X), nsamples)
    d <- dist(X[sam,])
    sigma <- 2 * min(d)
  }

  wfun <- if (weight_mode == "heat") {
    function(x) heat_kernel(x, sigma)
  } else if (weight_mode == "binary") {
    function(x) (x>0) * 1
  } else if (weight_mode == "normalized") {
    function(x) normalized_heat_kernel(x,sigma,ncol(X))
  }

  if (neighbor_mode == "knn") {
    W <- weighted_knn(X, k, FUN=wfun)
  } else if (neighbor_mode == "supervised") {
      assertthat::assert_that(!is.null(labels))
      labels <- as.factor(labels)
      if (k == 0) {
        W <- label_matrix2(labels, labels)
      } else if (k > 0) {

        M <- do.call(rbind, lapply(levels(labels), function(lev) {
          idx <- which(labels == lev)
          M <- weighted_knn(X[idx,], k, FUN=wfun)
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
      browser()
      idx <- which(labels != lev)
      M <- weighted_knn(X[idx,], k, FUN=wfun)
      Mind <- which(M > 0, arr.ind=TRUE)
      cbind(idx[Mind[,1]], idx[Mind[,2]], M[Mind])
    }))

    Matrix::sparseMatrix(i=M[,1], j=M[,2], x=M[,3], dims=c(nrow(X), nrow(X)), use.last.ij = TRUE)

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


#' weighted_knn
#'
#' @param X the data matrix, where rows are instances and columns are features
#' @param eps the e-radius used to threshold neighbor inclusion
#' @param FUN the function used to convert euclidean distances to similarities
#' @importFrom assertthat assert_that
#' @importFrom FNN get.knn
#' @importFrom Matrix t
#' @export
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

  nnd <- nn$nn.dist

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

