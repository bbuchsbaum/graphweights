#' @importFrom RSpectra eigs
#' @export
commute_time <- function(A, ncomp=nrow(A)-1) {
  D <- Matrix::rowSums(A)

  Dtilde <- Diagonal(x= D^(-1/2))

  M <- Dtilde %*% A %*% Dtilde

  decomp <- RSpectra::eigs_sym(M, k=ncomp+1)

  pii <- D/sum(A)
  v <- decomp$vectors[, 2:(ncomp+1)]
  ev <- decomp$values[2:(ncomp+1)]

  gap <- decomp$values[1] - decomp$values[2]

  cds <- sweep(v, 2, sqrt(1 - ev), "/")
  cds <- sweep(cds, 1, 1/sqrt(pii), "*")
  cds

  ret <- list(vectors=decomp$vectors, values=decomp$values, cds=cds, gap=gap)
  class(ret) <- c("commute_time", "list")
  ret
}

correct_commute <- function(A, cds) {
  #R = get_commute_distance_using_Laplacian(A)

  n = nrow(A)

  # compute commute time limit expression:
  d = rowSums(A)
  Rlimit <- t(matrix(1/d, nrow=n, ncol=n, byrow=TRUE)) + t(matrix(1/d, nrow=n, ncol=n))
  #Rlimit = np.tile((1. / d), (n, 1)).T + np.tile((1. / d), (n, 1))

  # compute correction term u_{ij}:
  tmp <- t(matrix(diag(A), nrow=n, ncol=n, byrow=TRUE)) / t(matrix(d, nrow=n, ncol=n))
  #tmp = np.tile(np.diag(A), (n, 1)).T / np.tile(d, (n, 1))
  tmp2 <- t(matrix(d, nrow=n, ncol=n, byrow=TRUE)) * matrix(d, nrow=n, ncol=n, byrow=TRUE)
  #tmp2 = np.tile(d, (n, 1)).T * np.tile(d, (n, 1))
  uAu = tmp + t(tmp) - 2 * A / tmp2

  # compute amplified commute:
  #tmp3 = R - Rlimit - uAu

  # enforce 0 diagonal:
  #D = tmp3 - np.diag(np.diag(tmp3))
}

commute_time_dist <- function(A, ncomp=nrow(A)-1) {
  n = nrow(A)
  D <- Matrix::Diagonal(x=rowSums(A))
  L = D - A

  dinv = Matrix::Diagonal(x=1/rowSums(A))

  #Linv = linalg.inv(L + np.ones(L.shape)/n) - np.ones(L.shape)/n
  Linv <- solve(L + matrix(1, nrow(A), ncol(A))/n)
  Linv_diag <- matrix(diag(Linv), n, 1)
  #Linv_diag = np.diag(Linv).reshape((n, 1))
  Rexact = Linv_diag %*% matrix(1, ncol=n) + matrix(1, nrow=n) %*% t(Linv_diag) - 2 * Linv
  #Rexact = Linv_diag * np.ones((1, n)) + np.ones((n, 1)) * Linv_diag.T - 2 * Linv
  # convert from a resistance distance to a commute time distance

  Rexact <- Rexact * sum(A)
}




# cds <- commute_time(A, ncomp=4)
# plot(cds$cds[,1], cds$cds[,2], type="p", col=iris$Species)
#
# out <- do.call(rbind, lapply(seq(4, 150, by=5), function(k) {
#   do.call(rbind, lapply(seq(.5, 2.5, by=.3), function(sigma) {
#     A=neighborweights::similarity_matrix(iris[,1:4],  k=k, sigma=sigma, weight_mode="heat")
#     ret <- commute_time(A, ncomp=4)
#     message("k = ", k, " : sigma = ", sigma, " : gap = ", ret$gap)
#     c(k=k, sigma=sigma, gap=ret$gap)
#   }))
# }))
