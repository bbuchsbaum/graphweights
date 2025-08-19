#' @useDynLib neighborweights, .registration=TRUE
#' @importFrom Rcpp sourceCpp
NULL

#' Compute a spatial autocorrelation matrix
#'
#' This function computes a spatial autocorrelation matrix using a radius-based nearest neighbor search.
#' The function leverages the mgcv package to fit a generalized additive model (GAM) to the data
#' and constructs the autocorrelation matrix using the fitted model.
#'
#' @param X A numeric matrix or data.frame, where each column represents a variable and each row represents an observation.
#' @param cds A numeric matrix or data.frame of spatial coordinates (x, y, or more dimensions) with the same number of rows as X.
#' @param radius A positive numeric value representing the search radius for the radius-based nearest neighbor search. Default is 8.
#' @param nsamples A positive integer indicating the number of samples to be taken for fitting the GAM. Default is 1000.
#'
#' @return A sparse matrix representing the spatial autocorrelation matrix for the input data.
#'
#' @importFrom mgcv gam
#' @importFrom Matrix Diagonal
#' @importFrom rflann RadiusSearch
#' @importFrom Matrix sparseMatrix
#' @importFrom assertthat assert_that
#' @export
spatial_autocor <- function(X, cds, radius=8, nsamples=1000) {
  assertthat::assert_that(radius > 0)
  if (nrow(cds) < nsamples) {
    nsamples <- nrow(cds)
  }

  nabe <- rflann::RadiusSearch(cds, cds,radius=radius)
  ids <- sample(1:nrow(nabe$indices), nsamples)

  ret <- do.call(rbind, lapply(ids, function(id) {
    ind <- nabe$indices[id,]
    d <- nabe$distances[id,]
    valmat <- X[,ind[-1]]
    vals <- X[,ind[1]]
    cvals <- cor(vals,valmat)
    data.frame(cor=as.vector(cvals), d=sqrt(d[-1]))
  }))

  gam.1 <- gam(cor ~ s(d), data=ret)

  trip <- do.call(rbind, lapply(1:nrow(cds), function(i) {
    i <- nabe$indices[i,1]
    j <- nabe$indices[i,-1]
    d <- sqrt(nabe$distances[i,-1])
    cbind(i,j,d)
  }))

  trip[,3] <- predict(gam.1, newdata=data.frame(d=trip[,3]))
  smat <- Matrix::sparseMatrix(i=trip[,1], j=trip[,2], x=trip[,3], dims=c(nrow(cds),
                                                                  nrow(cds)))

  smat
}


#' Compute a pairwise adjacency matrix for multiple graphs
#'
#' This function computes a pairwise adjacency matrix for multiple graphs with a given set of spatial coordinates
#' and feature vectors. The function takes two user-defined functions to compute within-graph and between-graph
#' similarity measures.
#'
#' @param Xcoords A list of numeric matrices or data.frames containing the spatial coordinates of the nodes of each graph.
#' @param Xfeats A list of numeric matrices or data.frames containing the feature vectors for the nodes of each graph.
#' @param fself A function that computes similarity for nodes within the same graph (e.g., Xi_1, Xi_2).
#' @param fbetween A function that computes similarity for nodes across graphs (e.g., Xi_1, Xj_1).
#'
#' @return A sparse matrix representing the pairwise adjacency matrix for the input graphs.
#'
#' @importFrom assertthat assert_that
#' @importFrom Matrix sparseMatrix
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

#' Compute the spatial Laplacian matrix of a coordinate matrix
#'
#' This function computes the spatial Laplacian matrix of a given coordinate matrix using specified parameters.
#'
#' @param coord_mat A numeric matrix representing coordinates
#' @param dthresh Numeric, the distance threshold for adjacency (default is 1.42)
#' @param nnk Integer, the number of nearest neighbors for adjacency (default is 27)
#' @param weight_mode Character, the mode for computing weights, either "binary" or "heat" (default is "binary")
#' @param sigma Numeric, the sigma parameter for the heat kernel (default is dthresh/2)
#' @param normalized Logical, whether the adjacency matrix should be normalized (default is TRUE)
#' @param stochastic Logical, whether the adjacency matrix should be stochastic (default is FALSE)
#'
#' @return A sparse symmetric matrix representing the computed spatial Laplacian
#'
#' @examples
#' # Create an example coordinate matrix
#' coord_mat <- matrix(c(1, 2, 3, 4, 5, 6), nrow = 3, ncol = 2)
#'
#' # Compute the spatial Laplacian matrix using the binary weight mode
#' result <- spatial_laplacian(coord_mat, dthresh = 1.42, nnk = 27, weight_mode = "binary")
#'
#' @export
spatial_laplacian <- function(coord_mat, dthresh=1.42, nnk=27,weight_mode=c("binary", "heat"),
                              sigma=dthresh/2,
                              normalized=TRUE, stochastic=FALSE) {

  weight_mode <- match.arg(weight_mode)
  adj <- spatial_adjacency(coord_mat, dthresh, nnk,weight_mode, sigma,
                           include_diagonal=FALSE, normalized, stochastic)
  L <- Diagonal(x=rowSums(adj))  - adj
}


#' Spatial Laplacian of Gaussian for coordinates
#'
#' This function computes the spatial Laplacian of Gaussian for a given coordinate matrix using a specified sigma value.
#'
#' @param coord_mat A numeric matrix representing coordinates
#' @param sigma Numeric, the sigma parameter for the Gaussian smoother (default is 2)
#'
#' @return A sparse symmetric matrix representing the computed spatial Laplacian of Gaussian
#'
#' @examples
#' # Create an example coordinate matrix
#' coord_mat <- matrix(c(1, 2, 3, 4, 5, 6), nrow = 3, ncol = 2)
#'
#' # Compute the spatial Laplacian of Gaussian with sigma = 2
#' result <- spatial_lap_of_gauss(coord_mat, sigma = 2)
#'
#' @export
spatial_lap_of_gauss <- function(coord_mat, sigma=2) {
  lap <- spatial_laplacian(coord_mat, weight_mode="binary", normalized=FALSE, nnk=ncol(coord_mat)^3, dthresh=sigma*2.5)
  #lap <- spatial_laplacian(coord_mat, weight_mode="binary", nnk=ncol(coord_mat)^3)
  adj <- spatial_smoother(coord_mat,  sigma=sigma)
  lap %*% adj
}

#' Compute the Difference of Gaussians for a coordinate matrix
#'
#' This function computes the Difference of Gaussians for a given coordinate matrix using specified sigma values and the number of nearest neighbors.
#'
#' @param coord_mat A numeric matrix representing coordinates
#' @param sigma1 Numeric, the first sigma parameter for the Gaussian smoother (default is 2)
#' @param sigma2 Numeric, the second sigma parameter for the Gaussian smoother (default is sigma1 * 1.6)
#' @param nnk Integer, the number of nearest neighbors for adjacency (default is 3^(ncol(coord_mat)))
#'
#' @return A sparse symmetric matrix representing the computed Difference of Gaussians
#'
#' @examples
#' # Create an example coordinate matrix
#' coord_mat <- matrix(c(1, 2, 3, 4, 5, 6), nrow = 3, ncol = 2)
#'
#' # Compute the Difference of Gaussians with sigma1 = 2 and sigma2 = 3.2
#' result <- difference_of_gauss(coord_mat, sigma1 = 2, sigma2 = 3.2)
#'
#' @export
difference_of_gauss <- function(coord_mat, sigma1=2, sigma2=sigma1 * 1.6, nnk=3^(ncol(coord_mat))) {
  adj1 <- spatial_smoother(coord_mat,  sigma=sigma1, stochastic=TRUE, nnk=nnk)
  adj2 <- spatial_smoother(coord_mat,  sigma=sigma2,stochastic=TRUE, nnk=nnk)
  adj1 - adj2
}





#' Compute the spatial smoother matrix for a coordinate matrix
#'
#' This function computes the spatial smoother matrix for a given coordinate matrix using specified parameters.
#'
#' @param coord_mat A numeric matrix representing coordinates
#' @param sigma Numeric, the sigma parameter for the Gaussian smoother (default is 5)
#' @param nnk Integer, the number of nearest neighbors for adjacency (default is 3^(ncol(coord_mat)))
#' @param stochastic Logical, whether the adjacency matrix should be doubly stochastic (default is TRUE)
#'
#' @return A sparse symmetric matrix representing the computed spatial smoother
#'
#' @examples
#' # Create an example coordinate matrix
#' coord_mat <- matrix(c(1, 2, 3, 4, 5, 6), nrow = 3, ncol = 2)
#'
#' # Compute the spatial smoother matrix with sigma = 5
#' result <- spatial_smoother(coord_mat, sigma = 5, nnk = 3^(ncol(coord_mat)), stochastic = TRUE)
#'
#' @export
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


#' Compute the spatial adjacency matrix for a coordinate matrix
#'
#' This function computes the spatial adjacency matrix for a given coordinate matrix using specified parameters.
#' Adjacency is determined by distance threshold and the maximum number of neighbors.
#'
#' @param coord_mat A numeric matrix representing the spatial coordinates
#' @param dthresh Numeric, the distance threshold defining the radius of the neighborhood (default is sigma*3)
#' @param nnk Integer, the maximum number of neighbors to include in each spatial neighborhood (default is 27)
#' @param weight_mode Character, the mode for computing weights, either "binary" or "heat" (default is "binary")
#' @param sigma Numeric, the bandwidth of the heat kernel if weight_mode == "heat" (default is 5)
#' @param include_diagonal Logical, whether to assign 1 to diagonal elements (default is TRUE)
#' @param normalized Logical, whether to make row elements sum to 1 (default is TRUE)
#' @param stochastic Logical, whether to make column elements also sum to 1 (only relevant if normalized == TRUE) (default is FALSE)
#'
#' @return A sparse symmetric matrix representing the computed spatial adjacency
#'
#' @importFrom rflann RadiusSearch
#' @importFrom Matrix sparseMatrix
#' @export
#'
#' @examples
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


#' Compute the doubly stochastic matrix from a given matrix
#'
#' This function iteratively computes the doubly stochastic matrix from a given input matrix.
#' A doubly stochastic matrix is a matrix in which both row and column elements sum to 1.
#'
#' @param A A numeric matrix for which to compute the doubly stochastic matrix
#' @param iter Integer, the number of iterations to perform (default is 30)
#'
#' @return A numeric matrix representing the computed doubly stochastic matrix
#'
#' @export
#'
#' @examples
#' # Create an example matrix
#' A <- matrix(c(2, 4, 6, 8, 10, 12), nrow = 3, ncol = 2)
#'
#' # Compute the doubly stochastic matrix
#' result <- make_doubly_stochastic(A, iter = 30)
make_doubly_stochastic <- function(A, iter=30) {
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

#' Normalize Adjacency Matrix
#'
#' This function normalizes an adjacency matrix by dividing each element by the product of the square root of the corresponding row and column sums.
#' Optionally, it can also symmetrize the normalized matrix by averaging it with its transpose.
#'
#' @param sm A sparse adjacency matrix representing the graph.
#' @param symmetric A logical value indicating whether to symmetrize the matrix after normalization (default: TRUE).
#'
#' @return A normalized and, if requested, symmetrized adjacency matrix.
#'
#' @examples
#' set.seed(123)
#' A <- matrix(runif(100), 10, 10)
#' A_normalized <- normalize_adjacency(A)
#'
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

#' Cross-adjacency matrix with feature weighting
#'
#' @param coord_mat1 the first coordinate matrix (the query)
#' @param coord_mat2 the second coordinate matrix (the reference)
#' @param feature_mat1 the first feature matrix
#' @param feature_mat2 the second feature matrix
#' @param wsigma the sigma for the feature heat kernel
#' @param nnk the maximum number of spatial nearest neighbors to include
#' @param maxk the maximum number of neighbors to include within spatial window
#' @param alpha the mixing weight for the spatial distance (1=all spatial weighting, 0=all feature weighting)
#' @param weight_mode the type of weighting to use: "binary" or "heat" (default: "binary")
#' @param sigma the spatial sigma for the heat kernel weighting (default: 1)
#' @param dthresh the threshold for the spatial distance
#' @param normalized whether to normalize the rows to sum to 1
#' @importFrom assertthat assert_that
#' @export
#'
#' @examples
#' set.seed(123)
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
                                    sigma, wsigma, alpha,
                                    maxk,
                                    weight_mode == "binary")

  sm <- sparseMatrix(i=triplet[,1], j=triplet[,2], x=triplet[,3],
                     dims=c(nrow(coord_mat1), nrow(coord_mat2)))

  if (normalized) {
    sm <- normalize_adjacency(sm)
  }

  sm

}



#' Bilateral Spatial Smoother
#'
#' This function computes a bilateral smoothing of the input data, which combines spatial and feature information to provide a smoothed representation of the data.
#'
#' @param coord_mat A matrix with the spatial coordinates of the data points, where each row represents a point and each column represents a coordinate dimension.
#' @param feature_mat A matrix with the feature vectors of the data points, where each row represents a point and each column represents a feature dimension.
#' @param nnk The number of nearest neighbors to consider for smoothing (default: 27).
#' @param s_sigma The spatial bandwidth in standard deviations (default: 2.5).
#' @param f_sigma The normalized feature bandwidth in standard deviations (default: 0.7).
#' @param stochastic A logical value indicating whether to make the resulting adjacency matrix doubly stochastic (default: FALSE).
#'
#' @return A sparse adjacency matrix representing the smoothed data.
#'
#' @examples
#' set.seed(123)
#' coord_mat <- as.matrix(expand.grid(1:10, 1:10))
#' feature_mat <- matrix(rnorm(100*10), 100, 10)
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
                              feature_mat, s_sigma, f_sigma)

  sm <- sparseMatrix(i=triplet[,1], j=triplet[,2], x=triplet[,3], dims=c(nrow(coord_mat), nrow(coord_mat)))

  if (stochastic) {
    sm <- make_doubly_stochastic(sm)
  }

  sm

}

#' Weighted Spatial Adjacency
#'
#' Constructs a spatial adjacency matrix, where weights are determined by a secondary feature matrix.
#'
#' @param coord_mat A matrix with the spatial coordinates of the data points, where each row represents a point and each column represents a coordinate dimension.
#' @param feature_mat A matrix with the feature vectors of the data points, where each row represents a point and each column represents a feature dimension. The number of rows in feature_mat must be equal to the number of rows in coord_mat.
#' @param wsigma The spatial weight scale (default: 0.73).
#' @param alpha The mixing parameter between 0 and 1 (default: 0.5). A value of 0 results in a purely spatial adjacency matrix, while a value of 1 results in a purely feature-based adjacency matrix.
#' @param nnk The number of nearest neighbors to consider (default: 27).
#' @param weight_mode The mode to use for weighting the adjacency matrix, either "binary" or "heat" (default: "binary").
#' @param sigma The bandwidth for heat kernel weights (default: 1).
#' @param dthresh The distance threshold for nearest neighbors (default: sigma * 2.5).
#' @param include_diagonal A logical value indicating whether to include diagonal elements in the adjacency matrix (default: TRUE).
#' @param normalized A logical value indicating whether to normalize the adjacency matrix (default: FALSE).
#' @param stochastic A logical value indicating whether to make the resulting adjacency matrix doubly stochastic (default: FALSE).
#'
#' @return A sparse adjacency matrix with weighted spatial relationships.
#'
#' @examples
#' \donttest{
#' set.seed(123)
#' coord_mat <- as.matrix(expand.grid(x=1:9, y=1:9, z=1:9))
#' fmat <- matrix(rnorm(nrow(coord_mat) * 100), nrow(coord_mat), 100)
#' wsa1 <- weighted_spatial_adjacency(coord_mat, fmat, nnk=3, weight_mode="binary", alpha=1, stochastic=TRUE)
#' wsa2 <- weighted_spatial_adjacency(coord_mat, fmat, nnk=27, weight_mode="heat", alpha=0, stochastic=TRUE, sigma=2.5)
#' }
#' @export
weighted_spatial_adjacency <- function(coord_mat, feature_mat, wsigma=.73, alpha=.5,
                                       nnk=27,
                                       weight_mode=c("binary", "heat"),
                                       sigma=1,
                                       dthresh=sigma*2.5,
                                       include_diagonal=TRUE,
                                       normalized=FALSE,
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
                              feature_mat, sigma, wsigma, alpha,
                              weight_mode == "binary")

  sm <- sparseMatrix(i=triplet[,1], j=triplet[,2], x=triplet[,3], dims=c(nrow(coord_mat), nrow(coord_mat)))

  if (!include_diagonal) {
    diag(sm) <- rep(0,nrow(sm))
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



#' Cross Spatial Adjacency
#'
#' Constructs a cross spatial adjacency matrix between two sets of points, where weights are determined by spatial relationships.
#'
#' @param coord_mat1 A matrix with the spatial coordinates of the first set of data points, where each row represents a point and each column represents a coordinate dimension.
#' @param coord_mat2 A matrix with the spatial coordinates of the second set of data points, where each row represents a point and each column represents a coordinate dimension.
#' @param dthresh The distance threshold for nearest neighbors (default: sigma * 3).
#' @param nnk The number of nearest neighbors to consider (default: 27).
#' @param weight_mode The mode to use for weighting the adjacency matrix, either "binary" or "heat" (default: "binary").
#' @param sigma The bandwidth for heat kernel weights (default: 5).
#' @param normalized A logical value indicating whether to normalize the adjacency matrix (default: TRUE).
#'
#' @return A sparse cross adjacency matrix with weighted spatial relationships between the two sets of points.
#'
#' @examples
#' \donttest{
#' coord_mat1 <- as.matrix(expand.grid(x=1:5, y=1:5, z=1:5))
#' coord_mat2 <- as.matrix(expand.grid(x=6:10, y=6:10, z=6:10))
#' csa <- cross_spatial_adjacency(coord_mat1, coord_mat2, nnk=3, weight_mode="binary", sigma=5, normalized=TRUE)
#' }
#'
#' @export
cross_spatial_adjacency <- function(coord_mat1, coord_mat2, dthresh=sigma*3,
                                    nnk=27, weight_mode=c("binary", "heat"),
                                    sigma=5,
                                    normalized=TRUE) {


  full_nn <- rflann::RadiusSearch(coord_mat1, coord_mat2, radius=dthresh^2,
                                  max_neighbour=nnk)

  weight_mode <- match.arg(weight_mode)

  nels <- sum(sapply(full_nn$indices, length))

  triplet <- spatial_weights(full_nn$indices, lapply(full_nn$distances, sqrt), sigma,
                             weight_mode == "binary")

  sm <- sparseMatrix(i=triplet[,1], j=triplet[,2], x=triplet[,3],
                     dims=c(nrow(coord_mat1), nrow(coord_mat2)))

  if (normalized) {
    D <- Diagonal(x=1/rowSums(sm))
    sm <- D %*% sm
  }

  sm
}

