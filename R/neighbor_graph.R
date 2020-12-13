
#' @export
neighbor_graph <- function(G, params=list(), type=NULL, classes=NULL, ...) {
  if (inherits(G, "Matrix") || inherits(G, "matrix")) {
    G <- igraph::graph_from_adjacency_matrix(G, mode="undirected", weighted=TRUE)
  }

  assert_that(inherits(G, "igraph"), msg="'G' must be one of tyhe types 'matrix', 'Matrix' or 'igraph' types")

  structure(list(
            G=G,
            params=params,
            ...),
            class=c(classes, "neighbor_graph"))
}


#' @export
adjacency.neighbor_graph <- function(x, attr="weight") {
  igraph::as_adjacency_matrix(W_g_mutual, attr=attr)
}

#' @export
laplacian.neighbor_graph <- function(x, normalized=FALSE) {
  A <- adjacency(x)
  D <- Matrix::Diagonal(x=rowSums(a))
  L <- D - A
  L
}

#' @export
neighbors.neighbor_graph <- function(x, i) {
  igraph::adjacent_vertices(x,i)
}

#' @export
node_density.neighbor_graph <- function(x, X) {
  L <- neighbors(x, 1:nrow(x))
  d <- sapply(1:length(L), function(i) {
    ind <- L[[i]]
    xi <- X[i,]
    xnn <- X[ind,,drop=FALSE]
    S <- sum(sweep(xnn, 2, xi, "-")^2)
    1/(length(ind)^2) * S
  })

  d
}

#' @export
non_neighbors.neighbor_graph <- function(x, i) {
  ret <- lapply(i, function(j) {
    which(x$G[j,] == 0)
  })
  names(ret) <- i
  ret
}

#' @export
nvertices.neighbor_graph <- function(x,...) {
  igraph::vcount(x$G)
}


