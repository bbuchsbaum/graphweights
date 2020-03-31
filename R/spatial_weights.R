
#' pairwise_adjacency
#'
#' @param Xcoords A list of coordinate matrices containing the spatial coordinates of the nodes of each graph
#' @param Xfeats A list of feature matrices containing the feature vectors for the nodes of each graph
#' @param fself a function that compute similarity for nodes of the same graph (e.g. Xi_1, Xi_2)
#' @param fbetween a function that compute similarity for nodes across graphs (e.g. Xi_1, Xj_1)
#' @export
pairwise_adjacency <- function(Xcoords, Xfeats, fself, fbetween) {
  assertthat::assert_that(length(Xcoords) == length(Xfeats))
  ninstances <- unlist(lapply(Xcoords, nrow))
  assertthat::assert_that(
    all(ninstances ==  unlist(lapply(Xfeats, nrow)))
  )

  assertthat::assert_that(is.function(fself))
  assertthat::assert_that(is.function(fbetween))

  nsets <- length(Xcoords)

  offsets <- cumsum(c(0, ninstances[1:(nsets-1)]))

  block_indices <- sapply(1:length(ninstances), function(i) {
    seq(offsets[i] + 1, offsets[i] + ninstances[i])
  })

  cmb <- t(combn(1:length(Xcoords), 2))

  bet <- do.call(rbind, lapply(1:nrow(cmb), function(i) {
    a <- cmb[i,1]
    b <- cmb[i,2]
    sm <- fbetween(Xcoords[[a]], Xcoords[[b]], Xfeats[[a]], Xfeats[[b]])
    sm_nc <- as (sm, "dgTMatrix")
    r1 <- cbind (i = sm_nc@i + 1 + offsets[a], j = sm_nc@j + 1 + offsets[b], x = sm_nc@x)
    r2 <- cbind (i = sm_nc@j + 1 + offsets[b], j = sm_nc@i + 1 + offsets[a], x = sm_nc@x)
    rbind(r1,r2)
  }))

  w <- do.call(rbind, lapply(1:length(Xcoords), function(i) {
    sm <- fself(Xcoords[[i]], Xfeats[[i]])
    sm_nc <- as (sm, "dgTMatrix")
    cbind (i = sm_nc@i + 1 + offsets[i], j = sm_nc@j + 1 + offsets[i], x = sm_nc@x)
  }))

  tot <- rbind(w,bet)
  ret <- sparseMatrix(i=tot[,1], j=tot[,2], x=tot[,3])
  ret

}

#' construct a laplacian matrix in a local spatial neighborhood
#'
#' @importFrom Matrix Diagonal
#' @export
#' @inheritParams spatial_adjacency
#' @examples
#'
#' coord_mat <- as.matrix(expand.grid(x=1:3, y=1:3, z=1:3))
#' lp <- spatial_laplacian(coord_mat)
#' all(dim(lp) == c(27,27))
spatial_laplacian <- function(coord_mat, dthresh=1.42, nnk=27,weight_mode=c("binary", "heat"), sigma=dthresh/2,
                              normalized=TRUE, stochastic=FALSE) {

  weight_mode <- match.arg(weight_mode)
  adj <- spatial_adjacency(coord_mat, dthresh, nnk,weight_mode, sigma,
                           include_diagonal=FALSE, normalized, stochastic)
  L <- Diagonal(x=rowSums(adj))  - adj
}


#' construct a laplacian-of-gaussian matrix in a local spatial neighborhood
#'
#' @param coord_mat the matrix of coordinates
#' @param sigma the standard deviation of the spatial smoother
#' @export
#' @examples
#'
#' coord_mat <- as.matrix(expand.grid(x=1:10, y=1:10))
#' lp <- spatial_lap_of_gauss(coord_mat)
#' all(dim(lp) == c(100,100))
spatial_lap_of_gauss <- function(coord_mat, sigma=2) {
  lap <- spatial_laplacian(coord_mat, weight_mode="binary", normalized=FALSE, nnk=ncol(coord_mat)^3, dthresh=sigma*2.5)
  #lap <- spatial_laplacian(coord_mat, weight_mode="binary", nnk=ncol(coord_mat)^3)
  adj <- spatial_smoother(coord_mat,  sigma=sigma)
  lap %*% adj
}


#' spatial_smoother
#'
#' @param coord_mat the matrix of coordinates
#' @param sigma the standard deviation of the smoother
#' @param nnk the number of neighbors in the kernel
#' @param stochastic make doubly stochastic
#' @export
#'
#' @examples
#' coord_mat <- as.matrix(expand.grid(x=1:10, y=1:10))
#' sm <- spatial_smoother(coord_mat, sigma=3)
#' sm2 <- spatial_smoother(coord_mat, sigma=6, stochastic=FALSE,nnk=50)
spatial_smoother <- function(coord_mat, sigma=5, nnk=3^(ncol(coord_mat)), stochastic=TRUE) {
  adj <- spatial_adjacency(coord_mat, dthresh=sigma*1, nnk=nnk,weight_mode="heat", sigma,
                           include_diagonal=TRUE, normalized=FALSE, stochastic=FALSE)

  D <- Matrix::rowSums(adj)

  adj <- 1/D * adj
  if (stochastic) {
    adj <- make_doubly_stochastic(adj)
    adj <- (adj + t(adj))/2
  }

  adj
}


#' spatial_adjacency
#'
#' @param coord_mat the matrix of spatial coordinates
#' @param dthresh the distance threshold defining the radius of the neighborhood
#' @param nnk the maximum number of neighbors to include in each spatial neighborhood
#' @param weight_mode a binary or heat kernel
#' @param sigma the bandwidth of the heat kernel if weight_mode == "heat"
#' @param include_diagonal diagonal elements assigned 1
#' @param normalized make row elements sum to 1
#' @param stochastic make column elements also sum to 1 (only relevant if \code{normalized == TRUE})
#' @importFrom rflann RadiusSearch
#' @importFrom Matrix sparseMatrix
#' @export
#'
#' @examples
#'
#' coord_mat = as.matrix(expand.grid(x=1:6, y=1:6))
#' sa <- spatial_adjacency(coord_mat)
spatial_adjacency <- function(coord_mat, dthresh=sigma*3, nnk=27, weight_mode=c("binary", "heat"), sigma=5,
                              include_diagonal=TRUE, normalized=TRUE, stochastic=FALSE) {

  sm <- cross_spatial_adjacency(coord_mat, coord_mat, dthresh=dthresh,
                                nnk=nnk, weight_mode=weight_mode,
                                sigma=sigma, normalized=FALSE)

  if (include_diagonal) {
    diag(sm) <- rep(1, nrow(sm))
  } else {
    diag(sm) <- rep(0, nrow(sm))
  }

  if (normalized) {
    sm <- normalize_adjacency(sm)
  }

  if (stochastic) {
    ## make doubly stochastic (row and columns sum to 1)
    sm <- make_doubly_stochastic(sm)
  }

  sm
}


#' make_doubly_stochastic
#'
#' @param A the smoothing matrix
#' @param iter number of iterations
#' @export
make_doubly_stochastic <- function(A, iter=15) {
  r <- rep(1, nrow(A))
  #Given A, let (n, n) = size(A) and initialize r = ones(n, 1);

  for (k in 1:iter) {
    c <- 1/t(A) %*% r
    r <- 1/(A %*% c)
  }

  C <- Diagonal(x=c[,1])
  R <- Diagonal(x=r[,1])

  Ab <- R %*% A %*% C
}

# make_doubly_stochastic <- function(A, iter=8) {
#   e <- 100000
#
#  for (i in 1:iter) {
#     D <- rowSums(A)
#     D1 <- 1/sqrt(D)
#     A <- D1 * A * D1
#     A <- (A + t(A))/2
#
#   }
#   A
#
# }

#' normalize an adjacency matrix
#'
#' @param sm the adjacency matrix
#' @param symmetric whether to symmetrize after normalizing
#' @importFrom Matrix rowSums
#' @export
normalize_adjacency <- function(sm, symmetric=TRUE) {
  D <- Matrix::rowSums(sm)

  D1 <- Diagonal(length(D), x=1/sqrt(D))
  sm <- D1 %*% sm %*% D1

  if (symmetric) {
    sm <- (sm + t(sm))/2
  }

  sm
}

#' Compute the adjacency between two adjacency matrices weighted by two corresponding feature matrices.
#'
#' @param coord_mat1 the first coordinate matrix (the query)
#' @param coord_mat2 the second coordinate matrix (the reference)
#' @param feature_mat1 the first feature matrix
#' @param feature_mat2 the second feature matrix
#' @param wsigma the sigma for the feature heat kernel
#' @param nnk the maximum number of spatial nearest neighbors to include
#' @param maxk the maximum number of neighbors to include within spatial window
#' @param alpha the mixing weight for the spatial distance (1=all spatial weighting, 0=all feature weighting)
#' @param dthresh the threshold for the spatial distance
#' @param normalized whether to normalize the rows to sum to 1
#' @importFrom assertthat assert_that
#' @export
#'
#' @examples
#'
#' coords <- as.matrix(expand.grid(1:5, 1:5))
#' fmat1 <- matrix(rnorm(5*25), 25, 5)
#' fmat2 <- matrix(rnorm(5*25), 25, 5)
#'
#' adj <- cross_weighted_spatial_adjacency(coords, coords, fmat1, fmat2)
cross_weighted_spatial_adjacency <- function(coord_mat1, coord_mat2,
                                             feature_mat1, feature_mat2,
                                             wsigma=.73, alpha=.5,
                                             nnk=27,
                                             maxk=nnk,
                                             weight_mode=c("binary", "heat"),
                                             sigma=1,
                                             dthresh=sigma*2.5,
                                             normalized=TRUE) {


  assert_that(nrow(feature_mat1) == nrow(coord_mat1))
  assert_that(nrow(feature_mat2) == nrow(coord_mat2))


  assert_that(alpha >= 0 && alpha <=1)
  assert_that(dthresh > 0)

  full_nn <- rflann::RadiusSearch(coord_mat1, coord_mat2, radius=dthresh^2,
                                  max_neighbour=nnk)

  weight_mode <- match.arg(weight_mode)

  lens <- sapply(full_nn$indices, length)
  lens <- ifelse(lens > maxk, maxk, lens)
  #nels <- sum(sapply(full_nn$indices, length))
  nels <- sum(lens)


  triplet <- cross_fspatial_weights(full_nn$indices, lapply(full_nn$distances, sqrt),
                                    feature_mat1, feature_mat2,
                                    nels, sigma, wsigma, alpha,
                                    maxk=maxk,
                                    weight_mode == "binary")

  sm <- sparseMatrix(i=triplet[,1], j=triplet[,2], x=triplet[,3],
                     dims=c(nrow(coord_mat1), nrow(coord_mat2)))

  if (normalized) {
    sm <- normalize_adjacency(sm)
  }

  sm

}



#' bilateral spatial smoother
#'
#' @param s_sigma spatial bandwidth in standard deviations
#' @param f_sigma normalized feature bandwidth in standard deviations
#' @inheritParams spatial_smoother
#'
#' @examples
#'
#' coord_mat = as.matrix(expand.grid(1:10, 1:10))
#' feature_mat = matrix(rnorm(100*10), 100, 10)
#'
#' S <- bilateral_smoother(coord_mat, feature_mat, nnk=8)
#'
#' @export
bilateral_smoother <- function(coord_mat, feature_mat, nnk=27, s_sigma=2.5, f_sigma=.7, stochastic=FALSE) {
  assert_that(nrow(feature_mat) == nrow(coord_mat))
  assert_that(nnk >= 4)

  ## find the set of spatial nearest neighbors
  full_nn <- rflann::RadiusSearch(coord_mat, coord_mat, radius=(s_sigma*2.5)^2, max_neighbour=nnk)

  nels <- sum(sapply(full_nn$indices, length))

  triplet <- bilateral_weights(full_nn$indices, lapply(full_nn$distances, sqrt),
                              feature_mat, nels, s_sigma, f_sigma)

  sm <- sparseMatrix(i=triplet[,1], j=triplet[,2], x=triplet[,3], dims=c(nrow(coord_mat), nrow(coord_mat)))

  if (stochastic) {
    sm <- make_doubly_stochastic(sm)
  }

  sm

}

#' construct a spatial adjacency matrix weighted by a secondary feature matrix
#'
#' @param coord_mat the coordinate matrix
#' @param feature_mat the feature matrix, where: (\code{nrow(feature_mat) == nrow(coord_mat)})
#' @inheritParams cross_weighted_spatial_adjacency
#' @importFrom assertthat assert_that
#' @export
#'
#' @examples
#'
#' coord_mat <- as.matrix(expand.grid(x=1:3, y=1:3))
#' fmat <- matrix(rnorm(9*3), 9, 3)
#' wsa1 <- weighted_spatial_adjacency(coord_mat, fmat, nnk=3, weight_mode="binary", alpha=1, stochastic=TRUE)
#' wsa2 <- weighted_spatial_adjacency(coord_mat, fmat, nnk=3, weight_mode="binary", alpha=0, stochastic=TRUE)
weighted_spatial_adjacency <- function(coord_mat, feature_mat, wsigma=.73, alpha=.5,
                                       nnk=27,
                                       weight_mode=c("binary", "heat"),
                                       sigma=1,
                                       dthresh=sigma*2.5,
                                       include_diagonal=TRUE,
                                       normalized=TRUE,
                                       stochastic=FALSE) {


  assert_that(nrow(feature_mat) == nrow(coord_mat))
  assert_that(alpha >= 0 && alpha <=1)
  assert_that(dthresh > 0)

  ## find the set of spatial nearest neighbors
  full_nn <- rflann::RadiusSearch(coord_mat, coord_mat, radius=dthresh^2, max_neighbour=nnk)
  weight_mode <- match.arg(weight_mode)
  alpha2 <- 1-alpha


  nels <- sum(sapply(full_nn$indices, length))

  triplet <- fspatial_weights(full_nn$indices, lapply(full_nn$distances, sqrt),
                              feature_mat, nels, sigma, wsigma, alpha,
                              weight_mode == "binary")

  sm <- sparseMatrix(i=triplet[,1], j=triplet[,2], x=triplet[,3], dims=c(nrow(coord_mat), nrow(coord_mat)))

  if (!include_diagonal) {
    diag(sm) <- rep(0,nrow(sm))
  }

  if (normalized) {
    sm <- normalize_adjacency(sm)
    if (stochastic) {
      ## make doubly stochastic (row and columns sum to 1)
      sm <- make_doubly_stochastic(sm)
    }
  }

  sm

}



#' cross_spatial_adjacency
#'
#' @param coord_mat1 the first coordinate matrix
#' @param coord_mat2 the second coordinate matrix
#' @inheritParams spatial_adjacency
#' @importFrom rflann RadiusSearch
#' @importFrom Matrix sparseMatrix
#' @export
cross_spatial_adjacency <- function(coord_mat1, coord_mat2, dthresh=sigma*3,
                                    nnk=27, weight_mode=c("binary", "heat"),
                                    sigma=5,
                                    normalized=TRUE) {


  full_nn <- rflann::RadiusSearch(coord_mat1, coord_mat2, radius=dthresh^2,
                                  max_neighbour=nnk)

  weight_mode <- match.arg(weight_mode)

  nels <- sum(sapply(full_nn$indices, length))

  triplet <- spatial_weights(full_nn$indices, lapply(full_nn$distances, sqrt), nels, sigma,
                             weight_mode == "binary")

  sm <- sparseMatrix(i=triplet[,1], j=triplet[,2], x=triplet[,3],
                     dims=c(nrow(coord_mat1), nrow(coord_mat2)))

  if (normalized) {
    D <- Diagonal(x=1/rowSums(sm))
    sm <- D %*% sm
  }

  sm
}

