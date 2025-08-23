
#' Compute Markov diffusion kernel via eigen decomposition
#'
#' Efficient computation of the Markov diffusion kernel for a graph represented by
#' a sparse adjacency matrix. For large graphs, uses RSpectra to compute only the
#' leading k eigenpairs of the normalized transition matrix.
#'
#' @param A Square sparse adjacency matrix (dgCMatrix) of an undirected, weighted graph with non-negative entries.
#' @param t Diffusion time parameter (positive scalar).
#' @param k Number of leading eigenpairs to compute. If NULL, performs full eigendecomposition.
#' @param symmetric If TRUE (default), uses symmetric normalization to guarantee real eigenvalues.
#' @return dgCMatrix representing the diffusion kernel matrix.
#' @examples
#' library(Matrix)
#' # Create a simple 5x5 adjacency matrix (path graph)
#' A <- sparseMatrix(i = c(1, 2, 3, 4), j = c(2, 3, 4, 5), 
#'                   x = c(1, 1, 1, 1), dims = c(5, 5))
#' A <- A + t(A)  # Make symmetric
#' 
#' # Compute diffusion kernel with t = 0.5
#' K <- compute_diffusion_kernel(A, t = 0.5)
#' 
#' # Use only top 3 eigenpairs for efficiency
#' K_approx <- compute_diffusion_kernel(A, t = 0.5, k = 3)
#' @importFrom Matrix Diagonal crossprod
#' @importFrom RSpectra eigs
#' @export
compute_diffusion_kernel <- function(A, t, k = NULL, symmetric = TRUE) {
  if (!inherits(A, "dgCMatrix")) A <- as(A, "CsparseMatrix")
  n <- nrow(A)
  if (ncol(A) != n) stop("A must be square.")

  # Degree vector
  d <- Matrix::rowSums(A)
  if (any(d == 0)) warning("Isolated nodes present; they will remain isolated in kernel.")

  # Build transition matrix P
  if (symmetric) {
    Dm12 <- Diagonal(x = 1 / sqrt(d))
    P <- Dm12 %*% A %*% Dm12
  } else {
    Dinv <- Diagonal(x = 1 / d)
    P <- Dinv %*% A
  }

  # Eigen decomposition
  if (!is.null(k) && k < n) {
    # compute top k eigenpairs
    eig <- RSpectra::eigs(P, k = k, which = "LM")
    U <- eig$vectors    # n×k
    L <- eig$values     # length k
  } else {
    eig <- eigen(as.matrix(P), symmetric = symmetric)
    U <- eig$vectors    # n×n
    L <- eig$values     # length n
  }

  # Diffusion kernel: K = U diag(L^t) U^T
  Lt <- L^t
  K <- U %*% Diagonal(x = Lt) %*% t(U)

  # Return sparse
  return(as(K, "CsparseMatrix"))
}

#' Diffusion map embedding and distance
#'
#' Computes the diffusion map embedding of a graph and the pairwise diffusion distances
#' based on the leading eigenvectors of the normalized transition matrix.
#'
#' @param A Square sparse adjacency matrix (dgCMatrix) of an undirected, weighted graph.
#' @param t Diffusion time parameter (positive scalar).
#' @param k Number of diffusion coordinates to compute, excluding the trivial first coordinate.
#' @return A list with two components:
#'   \item{embedding}{n×k matrix of diffusion coordinates where n is the number of nodes.}
#'   \item{distances}{n×n matrix of squared diffusion distances between all node pairs.}
#' @examples
#' library(Matrix)
#' # Create a simple 4x4 adjacency matrix (path graph)
#' A <- sparseMatrix(i = c(1, 2, 3), j = c(2, 3, 4), x = c(1, 1, 1), dims = c(4, 4))
#' A <- A + t(A)  # Make symmetric
#' 
#' # Compute diffusion map with 2 coordinates
#' result <- compute_diffusion_map(A, t = 1.0, k = 2)
#' 
#' # View the embedding coordinates
#' print(result$embedding)
#' 
#' # Check diffusion distances
#' print(result$distances[1, ])  # distances from node 1
#' @importFrom RSpectra eigs
#' @importFrom Matrix Diagonal
#' @export
compute_diffusion_map <- function(A, t, k = 10) {
  if (!inherits(A, "dgCMatrix")) A <- as(A, "CsparseMatrix")
  n <- nrow(A)
  if (ncol(A) != n) stop("A must be square.")

  # Build symmetric transition P
  d <- Matrix::rowSums(A)
  Dm12 <- Diagonal(x = 1 / sqrt(d))
  P <- Dm12 %*% A %*% Dm12

  # Compute k+1 eigenpairs (first eigenvalue = 1)
  eig <- RSpectra::eigs(P, k = k + 1, which = "LM")
  lambdas <- eig$values
  U <- eig$vectors  # n × (k+1)

  # Sort in descending order
  ord <- order(Re(lambdas), decreasing = TRUE)
  lambdas <- Re(lambdas[ord])[2:(k+1)]
  U <- Re(U[, ord])[ ,2:(k+1), drop=FALSE]

  # Diffusion coordinates: psi_j = lambda_j^t * U[,j]
  coords <- t( t(U) * (lambdas^t) )

  # Squared diffusion distances: D_ij = sum_j (psi_j(i) - psi_j(j))^2
  # We can compute as: rowSums((coords[i,]-coords[j,])^2)
  # Efficient via crossprod: use matrix operations
  D2 <- as.matrix(dist(coords))^2

  list(embedding = coords, distances = D2)
}
