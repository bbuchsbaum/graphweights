

neighbord_graph <- function(G, neighbor_indices, params=list(), type=NULL, classes=NULL, ...) {
  structure(list(
            G=G,
            neighbor_indices=neighbor_indices,
            params=parms,
            type="weights",
            ...),
            class=c(classes, "weight_matrix"))
}

laplacian.weight_matrix <- function(x, normalized=FALSE) {
  D <- Matrix::Diagonal(x=rowSums(x$G))
  L <- x$D - x$G
  L
}

neighbors.weight_matrix <- function(x, i) {
  x$neighbor_indices[i]
}

neighbor_weights.weight_matrix <- function(x, i) {
  lapply(i, function(j) {
    ind <- x$neighbor_indices[[j]]
    G[ind,]
  })
}

non_neighbors <- function(x, i) {
  ind <- 1:ncol(x$G)
  lapply(i, function(j) {
    ind[-x$neighbor_indices[[j]]]
  })
}

nrow.weight_matrix <- function(x) nrow(x$W)
ncol.weight_matrix <- function(x) nrow(x$G)

