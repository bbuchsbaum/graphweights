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
