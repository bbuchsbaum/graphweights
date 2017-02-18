

#' @importFrom Matrix Diagonal
#' @export
spatial_laplacian <- function(coord_mat, dthresh=1.42, nnk=27) {
  adj <- spatial_adjacency(coord_mat=coord_mat, dthresh=dthresh, nnk=nnk)
  Diagonal(x=rowSums(adj))  - adj
}


#' spatial_adjacency
#'
#' @importFrom FNN get.knnx
#' @importFrom Matrix sparseMatrix
#' @export
spatial_adjacency <- function(coord_mat, dthresh=1.42, nnk=27) {
  full_nn <- FNN::get.knnx(coord_mat, coord_mat, nnk)

  triplet <- do.call(rbind, lapply(1:nrow(coord_mat), function(i) {
    dvals <- full_nn$nn.dist[i,1:nnk]
    keep <- dvals > 0 & dvals < dthresh
    if (keep[2]) {
      cbind(i,full_nn$nn.index[i,keep], 1)
    } else {
      NULL
    }
  }))

  sparseMatrix(i=triplet[,1], j=triplet[,2], x=triplet[,3])

}
