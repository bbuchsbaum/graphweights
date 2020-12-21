
#' @export
neighbor_graph.igraph <- function(x, params=list(), type=NULL, classes=NULL, ...) {

  structure(list(
            G=x,
            params=params,
            ...),
            class=c(classes, "neighbor_graph"))
}

#' @export
neighbor_graph.Matrix <- function(x, params=list(), type=NULL, classes=NULL, ...) {

  G <- igraph::graph_from_adjacency_matrix(x, mode="undirected", weighted=TRUE)

  structure(list(
    G=G,
    params=params,
    ...),
    class=c(classes, "neighbor_graph"))
}


#' @export
adjacency.neighbor_graph <- function(x, attr="weight") {
  igraph::as_adjacency_matrix(x$G, attr=attr)
}

#' @export
laplacian.neighbor_graph <- function(x, normalized=FALSE) {
  ## TODO normalize if requested
  A <- adjacency(x)
  D <- Matrix::Diagonal(x=rowSums(a))
  L <- D - A
  L
}

#' @export
neighbors.neighbor_graph <- function(x, i) {
  if (!missing(i)) {
    igraph::adjacent_vertices(x$G,i)
  } else {
    igraph::adjacent_vertices(x$G)
  }

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
edges.neighbor_graph <- function(x) {
  ret <- igraph::as_data_frame(cg$G, what="edges")
  cbind(ret[,1], ret[,2])
}


non_neighbors.neighbor_graph <- function(x, i) {
  #ret <- lapply(i, function(j) {
  #  which(x$G[j,] == 0)
  #})
  #names(ret) <- i
  #ret

  stop("not implemented")
}

#' @export
nvertices.neighbor_graph <- function(x,...) {
  igraph::vcount(x$G)
}


