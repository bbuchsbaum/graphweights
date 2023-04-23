
#' Neighbor Graph for igraph Objects
#'
#' Create a neighbor_graph object from an igraph object.
#'
#' @param x An igraph object.
#' @param params A list of parameters.
#' @param type A character string specifying the type of graph.
#' @param classes A character vector specifying additional classes for the object.
#' @param ... Additional arguments.
#'
#' @return A neighbor_graph object.
#'
#' @export
neighbor_graph.igraph <- function(x, params=list(), type=NULL, classes=NULL, ...) {

  structure(list(
            G=x,
            params=params,
            ...),
            class=c(classes, "neighbor_graph"))
}


#' Neighbor Graph for Matrix Objects
#'
#' Create a neighbor_graph object from a Matrix object.
#'
#' @param x A Matrix object.
#' @param params A list of parameters.
#' @param type A character string specifying the type of graph.
#' @param classes A character vector specifying additional classes for the object.
#' @param ... Additional arguments.
#'
#' @return A neighbor_graph object.
#'
#' @export
neighbor_graph.Matrix <- function(x, params=list(), type=NULL, classes=NULL, ...) {

  G <- igraph::graph_from_adjacency_matrix(x, mode="undirected", weighted=TRUE)

  structure(list(
    G=G,
    params=params,
    ...),
    class=c(classes, "neighbor_graph"))
}


#' Adjacency Matrix for neighbor_graph Objects
#'
#' Extract the adjacency matrix from a neighbor_graph object.
#'
#' @param x A neighbor_graph object.
#' @param attr A character string specifying the attribute to use for the adjacency matrix.
#'
#' @return A Matrix object representing the adjacency matrix.
#'
#' @export
adjacency.neighbor_graph <- function(x, attr="weight") {
  igraph::as_adjacency_matrix(x$G, attr=attr)
}



#' Laplacian Matrix for neighbor_graph Objects
#'
#' Compute the Laplacian matrix of a neighbor_graph object.
#'
#' @param x A neighbor_graph object.
#' @param normalized A logical value indicating whether the normalized Laplacian should be computed.
#'
#' @return A Matrix object representing the Laplacian matrix.
#'
#' @export
laplacian.neighbor_graph <- function(x, normalized=FALSE) {
  ## TODO normalize if requested
  A <- adjacency(x)
  D <- Matrix::Diagonal(x=rowSums(a))
  L <- D - A
  L
}

#' Neighbors for neighbor_graph Objects
#'
#' Retrieve the neighbors of a specific node or all nodes in a neighbor_graph object.
#'
#' @param x A neighbor_graph object.
#' @param i An integer specifying the index of the node for which neighbors should be retrieved.
#'
#' @return A list of neighbors.
#'
#' @export
neighbors.neighbor_graph <- function(x, i) {
  if (!missing(i)) {
    igraph::adjacent_vertices(x$G,i)
  } else {
    igraph::adjacent_vertices(x$G)
  }

}

#' Node Density for neighbor_graph Objects
#'
#' Compute the node density for a neighbor_graph object.
#'
#' @param x A neighbor_graph object.
#' @param X A data matrix containing the data points.
#'
#' @return A numeric vector containing the node densities.
#'
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

#' Edges for neighbor_graph Objects
#'
#' Retrieve the edges of a neighbor_graph object.
#'
#' @param x A neighbor_graph object.
#'
#' @return A matrix containing the edges of the neighbor_graph object.
#'
#' @export
edges.neighbor_graph <- function(x) {
  ret <- igraph::as_data_frame(cg$G, what="edges")
  cbind(ret[,1], ret[,2])
}


#' Non-neighbors for neighbor_graph Objects
#'
#' Retrieve the non-neighboring nodes of a given node in a neighbor_graph object.
#'
#' @param x A neighbor_graph object.
#' @param i The index of the node for which non-neighboring nodes will be returned.
#'
#' @return A list of node indices that are not neighbors of the given node.
#'
#' @export
non_neighbors.neighbor_graph <- function(x, i) {
  # Get all the vertices in the graph
  all_vertices <- seq_len(igraph::vcount(x$G))

  # Get the neighbors of the node 'i'
  neighbors_i <- neighbors(x, i)

  # Find the non-neighbors by excluding the neighbors and the node itself
  non_neighbors <- setdiff(all_vertices, c(i, neighbors_i))

  non_neighbors
}

#' @export
nvertices.neighbor_graph <- function(x,...) {
  igraph::vcount(x$G)
}


