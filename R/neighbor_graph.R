#' Create neighbor_graph from igraph or Matrix objects
#'
#' Create a neighbor_graph object from an igraph object or Matrix object.
#'
#' @param x An igraph or Matrix object.
#' @param params A list of parameters (default: empty list).
#' @param type A character string specifying the type of graph (currently unused).
#' @param classes A character vector specifying additional classes for the object.
#' @param ... Additional arguments.
#'
#' @return A neighbor_graph object wrapping the input graph structure.
#'
#' @examples
#' # Create from igraph
#' library(igraph)
#' g <- make_ring(5)
#' ng1 <- neighbor_graph(g)
#'
#' # Create from Matrix
#' adj_matrix <- Matrix::Matrix(c(0, 1, 1, 0, 1, 0, 1, 0, 0), 
#'                             nrow = 3, byrow = TRUE, sparse = TRUE)
#' ng2 <- neighbor_graph(adj_matrix)
#'
#' @rdname neighbor_graph
#' @method neighbor_graph igraph
#' @export
neighbor_graph.igraph <- function(x, params=list(), type=NULL, classes=NULL, ...) {

  structure(list(
            G=x,
            params=params,
            ...),
            class=c(classes, "neighbor_graph"))
}


#' @rdname neighbor_graph
#' @method neighbor_graph Matrix
#' @export
neighbor_graph.Matrix <- function(x, params=list(), type=NULL, classes=NULL, ...) {

  G <- igraph::graph_from_adjacency_matrix(x, mode="undirected", weighted=TRUE)

  structure(list(
    G=G,
    params=params,
    ...),
    class=c(classes, "neighbor_graph"))
}


#' Extract adjacency matrix from neighbor_graph object
#'
#' Extract the adjacency matrix from a neighbor_graph object.
#'
#' @param x A neighbor_graph object.
#' @param attr A character string specifying the edge attribute to use for weights (default: "weight").
#' @param ... Additional arguments (currently ignored).
#'
#' @return A sparse Matrix object representing the adjacency matrix.
#'
#' @examples
#' # Create example neighbor_graph
#' adj_matrix <- Matrix::Matrix(c(0, 1, 1, 0, 1, 0, 1, 0, 0), 
#'                             nrow = 3, byrow = TRUE, sparse = TRUE)
#' ng <- neighbor_graph(adj_matrix)
#' adj <- adjacency(ng)
#'
#' @method adjacency neighbor_graph
#' @export
adjacency.neighbor_graph <- function(x, attr="weight", ...) {
  igraph::as_adjacency_matrix(x$G, attr=attr)
}



#' Compute Laplacian matrix for neighbor_graph object
#'
#' Compute the Laplacian matrix of a neighbor_graph object.
#'
#' @param x A neighbor_graph object.
#' @param normalized A logical value indicating whether the normalized Laplacian should be computed (default: FALSE).
#' @param ... Additional arguments (currently ignored).
#'
#' @return A sparse Matrix object representing the Laplacian matrix.
#'
#' @details
#' The unnormalized Laplacian is computed as L = D - A where D is the degree matrix
#' and A is the adjacency matrix. The normalized Laplacian is computed as 
#' L_sym = I - D^(-1/2) A D^(-1/2).
#'
#' @examples
#' adj_matrix <- Matrix::Matrix(c(0, 1, 1, 0, 1, 0, 1, 0, 0), 
#'                             nrow = 3, byrow = TRUE, sparse = TRUE)
#' ng <- neighbor_graph(adj_matrix)
#' L <- laplacian(ng)
#' L_norm <- laplacian(ng, normalized = TRUE)
#'
#' @method laplacian neighbor_graph
#' @export
laplacian.neighbor_graph <- function(x, normalized=FALSE, ...) {
  A <- adjacency(x)
  N <- nrow(A)
  
  if (normalized) {
    d <- Matrix::rowSums(A)
    # Add epsilon to prevent division by zero for isolated nodes
    d_inv_sqrt <- ifelse(d > 0, 1 / sqrt(d), 0)
    D_inv_sqrt <- Matrix::Diagonal(x = d_inv_sqrt)
    # Calculate I - D^{-1/2} A D^{-1/2}
    L_sym <- Matrix::Diagonal(N) - D_inv_sqrt %*% A %*% D_inv_sqrt
    # Ensure symmetry numerically
    L_sym <- (L_sym + Matrix::t(L_sym)) / 2
    return(L_sym)
  } else {
    # Compute unnormalized Laplacian L = D - A
    d <- Matrix::rowSums(A) # Correct variable name here
    D <- Matrix::Diagonal(x=d)
    L <- D - A
    return(L)
  }
}

#' Get neighbors for neighbor_graph object
#'
#' Retrieve the neighbors of a specific node or all nodes in a neighbor_graph object.
#'
#' @param x A neighbor_graph object.
#' @param i An integer specifying the index of the node for which neighbors should be retrieved. If missing, returns neighbors for all nodes.
#' @param ... Additional arguments (currently ignored).
#'
#' @return If i is provided, a list containing the neighbors of node i. If i is missing, a list with neighbors for all nodes.
#'
#' @examples
#' adj_matrix <- Matrix::Matrix(c(0, 1, 1, 0, 1, 0, 1, 0, 0), 
#'                             nrow = 3, byrow = TRUE, sparse = TRUE)
#' ng <- neighbor_graph(adj_matrix)
#' neighbors(ng, 1)  # Neighbors of node 1
#' neighbors(ng)     # All neighbors
#'
#' @method neighbors neighbor_graph
#' @export
neighbors.neighbor_graph <- function(x, i, ...) {
  if (!missing(i)) {
    igraph::adjacent_vertices(x$G,i)
  } else {
    igraph::adjacent_vertices(x$G)
  }

}

#' Compute node density for neighbor_graph object
#'
#' Compute the node density for a neighbor_graph object based on local neighborhoods.
#'
#' @param x A neighbor_graph object.
#' @param X A data matrix containing the data points, with rows as observations.
#' @param ... Additional arguments (currently ignored).
#'
#' @return A numeric vector containing the node densities.
#'
#' @details
#' Node density is computed as the average squared distance from each node to its neighbors,
#' normalized by the square of the neighborhood size.
#'
#' @examples
#' \donttest{
#' # Create example data and graph
#' X <- matrix(rnorm(30), nrow = 10, ncol = 3)
#' adj_matrix <- Matrix::Matrix(diag(10), sparse = TRUE)
#' ng <- neighbor_graph(adj_matrix)
#' densities <- node_density(ng, X)
#' }
#'
#' @method node_density neighbor_graph
#' @export
node_density.neighbor_graph <- function(x, X, ...) {
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

#' Extract edges from neighbor_graph object
#'
#' Retrieve the edges of a neighbor_graph object.
#'
#' @param x A neighbor_graph object.
#' @param ... Additional arguments (currently ignored).
#'
#' @return A two-column matrix containing the edges, where each row represents an edge between two nodes.
#'
#' @examples
#' adj_matrix <- Matrix::Matrix(c(0, 1, 1, 0, 1, 0, 1, 0, 0), 
#'                             nrow = 3, byrow = TRUE, sparse = TRUE)
#' ng <- neighbor_graph(adj_matrix)
#' edge_list <- edges(ng)
#'
#' @method edges neighbor_graph
#' @export
edges.neighbor_graph <- function(x, ...) {
  ret <- igraph::as_data_frame(x$G, what="edges")
  cbind(ret[,1], ret[,2])
}


#' Get non-neighbors for neighbor_graph object
#'
#' Retrieve the non-neighboring nodes of a given node in a neighbor_graph object.
#'
#' @param x A neighbor_graph object.
#' @param i The index of the node for which non-neighboring nodes will be returned.
#' @param ... Additional arguments (currently ignored).
#'
#' @return A numeric vector of node indices that are not neighbors of the given node (excluding the node itself).
#'
#' @examples
#' adj_matrix <- Matrix::Matrix(c(0, 1, 1, 0, 1, 0, 1, 0, 0), 
#'                             nrow = 3, byrow = TRUE, sparse = TRUE)
#' ng <- neighbor_graph(adj_matrix)
#' non_neighbors(ng, 1)  # Non-neighbors of node 1
#'
#' @method non_neighbors neighbor_graph
#' @export
non_neighbors.neighbor_graph <- function(x, i, ...) {
  # Get all the vertices in the graph
  all_vertices <- seq_len(igraph::vcount(x$G))

  # Get the neighbors of the node 'i'
  neighbors_i <- neighbors(x, i)

  # Find the non-neighbors by excluding the neighbors and the node itself
  non_neighbors <- setdiff(all_vertices, c(i, neighbors_i))

  non_neighbors
}

#' Get number of vertices in neighbor_graph object
#'
#' Get the number of vertices in a neighbor_graph object.
#'
#' @param x A neighbor_graph object.
#' @param ... Additional arguments (currently ignored).
#'
#' @return An integer representing the number of vertices in the graph.
#'
#' @examples
#' adj_matrix <- Matrix::Matrix(c(0, 1, 1, 0, 1, 0, 1, 0, 0), 
#'                             nrow = 3, byrow = TRUE, sparse = TRUE)
#' ng <- neighbor_graph(adj_matrix)
#' nvertices(ng)  # Should return 3
#'
#' @method nvertices neighbor_graph
#' @export
nvertices.neighbor_graph <- function(x,...) {
  igraph::vcount(x$G)
}


