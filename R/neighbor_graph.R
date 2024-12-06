#' Neighbor Graph for igraph Objects
#'
#' @description
#' Creates a neighbor_graph object from an existing igraph object, providing a consistent
#' interface for graph-based operations in the graphweights package.
#'
#' @details
#' The neighbor_graph class provides a unified interface for working with neighborhood
#' relationships in graph structures. This constructor creates a neighbor_graph from
#' an existing igraph object, preserving its structure while adding functionality
#' specific to neighborhood-based operations.
#'
#' @param x An igraph object representing the graph structure
#' @param params A list of parameters used in graph construction (default: empty list)
#' @param type A character string specifying the type of graph (optional)
#' @param classes A character vector specifying additional classes to be added to
#'   the object (optional)
#' @param ... Additional arguments to be stored in the neighbor_graph object
#'
#' @return A neighbor_graph object with the following components:
#'   * G: The underlying igraph object
#'   * params: List of parameters used in graph construction
#'   * Additional components specified via ...
#'
#' @examples
#' \dontrun{
#' # Create a simple igraph object
#' library(igraph)
#' g <- make_ring(10)
#' 
#' # Convert to neighbor_graph
#' ng <- neighbor_graph(g)
#' }
#'
#' @seealso
#' * \code{\link{neighbor_graph.Matrix}} for creating from adjacency matrices
#' * \code{\link{adjacency.neighbor_graph}} for extracting adjacency information
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
#' @description
#' Creates a neighbor_graph object from a Matrix object representing an adjacency matrix.
#'
#' @details
#' Converts an adjacency matrix into a neighbor_graph object by first creating an
#' intermediate igraph representation. The matrix is interpreted as an undirected,
#' weighted graph where non-zero entries indicate edges and their weights.
#'
#' @param x A Matrix object representing an adjacency matrix
#' @param params A list of parameters used in graph construction (default: empty list)
#' @param type A character string specifying the type of graph (optional)
#' @param classes A character vector specifying additional classes to be added to
#'   the object (optional)
#' @param ... Additional arguments to be stored in the neighbor_graph object
#'
#' @return A neighbor_graph object with the following components:
#'   * G: An igraph object created from the adjacency matrix
#'   * params: List of parameters used in graph construction
#'   * Additional components specified via ...
#'
#' @examples
#' # Create a simple adjacency matrix
#' adj_mat <- Matrix::sparseMatrix(
#'   i = c(1,2,2,3),
#'   j = c(2,1,3,2),
#'   x = c(1,1,1,1),
#'   dims = c(3,3)
#' )
#' 
#' # Convert to neighbor_graph
#' ng <- neighbor_graph(adj_mat)
#'
#' @seealso
#' * \code{\link{neighbor_graph.igraph}} for creating from igraph objects
#' * \code{\link{adjacency.neighbor_graph}} for extracting adjacency information
#'
#' @importFrom igraph graph_from_adjacency_matrix
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
#' @description
#' Extracts the adjacency matrix representation from a neighbor_graph object,
#' optionally using a specific edge attribute for weights.
#'
#' @details
#' The adjacency matrix represents the graph structure where each entry (i,j)
#' indicates the presence and weight of an edge between vertices i and j.
#' For unweighted graphs, non-zero entries are 1. For weighted graphs,
#' the entries correspond to the specified edge attribute.
#'
#' @param x A neighbor_graph object
#' @param attr A character string specifying the edge attribute to use as weights
#'   (default: "weight")
#'
#' @return A sparse Matrix object representing the adjacency matrix where:
#'   * Non-zero entries indicate edges
#'   * Values correspond to edge weights (if weighted)
#'   * The matrix is symmetric for undirected graphs
#'
#' @examples
#' \dontrun{
#' # Create a neighbor graph
#' g <- make_ring(5)
#' ng <- neighbor_graph(g)
#' 
#' # Extract adjacency matrix
#' adj_mat <- adjacency(ng)
#' }
#'
#' @seealso
#' * \code{\link{laplacian.neighbor_graph}} for computing the Laplacian matrix
#' * \code{\link{neighbors.neighbor_graph}} for finding neighboring vertices
#'
#' @export
adjacency.neighbor_graph <- function(x, attr="weight") {
  igraph::as_adjacency_matrix(x$G, attr=attr)
}



#' Laplacian Matrix for neighbor_graph Objects
#'
#' @description
#' Computes the Laplacian matrix of a neighbor_graph object, with optional
#' normalization.
#'
#' @details
#' The Laplacian matrix L is computed as L = D - A, where:
#' * D is the degree matrix (diagonal matrix of vertex degrees)
#' * A is the adjacency matrix
#'
#' When normalized=TRUE, the normalized Laplacian L_norm = D^(-1/2) L D^(-1/2)
#' is computed instead.
#'
#' The Laplacian matrix is fundamental in spectral graph theory and is used for:
#' * Spectral clustering
#' * Graph partitioning
#' * Dimensionality reduction
#' * Computing graph properties
#'
#' @param x A neighbor_graph object
#' @param normalized Logical; whether to compute the normalized Laplacian
#'   (default: FALSE)
#'
#' @return A sparse Matrix object representing the Laplacian matrix
#'
#' @examples
#' \dontrun{
#' # Create a neighbor graph
#' g <- make_ring(5)
#' ng <- neighbor_graph(g)
#' 
#' # Compute standard Laplacian
#' L <- laplacian(ng)
#' 
#' # Compute normalized Laplacian
#' L_norm <- laplacian(ng, normalized = TRUE)
#' }
#'
#' @seealso
#' * \code{\link{adjacency.neighbor_graph}} for the adjacency matrix
#'
#' @importFrom Matrix Diagonal
#' @export
laplacian.neighbor_graph <- function(x, normalized=FALSE) {
  ## TODO normalize if requested
  A <- adjacency(x)
  D <- Matrix::Diagonal(x=rowSums(A))
  L <- D - A
  
  if (normalized) {
    D_inv_sqrt <- Matrix::Diagonal(x=1/sqrt(rowSums(A)))
    L <- D_inv_sqrt %*% L %*% D_inv_sqrt
  }
  
  L
}

#' Neighbors for neighbor_graph Objects
#'
#' @description
#' Retrieves the neighbors of a specific vertex or all vertices in a neighbor_graph
#' object.
#'
#' @details
#' For a given vertex i, returns all vertices that share an edge with i.
#' If no vertex is specified, returns a list of neighbors for all vertices
#' in the graph.
#'
#' @param x A neighbor_graph object
#' @param i Optional; integer specifying the vertex index for which to find neighbors.
#'   If missing, returns neighbors for all vertices
#'
#' @return If i is specified:
#'   * A vector of vertex indices representing neighbors of vertex i
#' If i is missing:
#'   * A list where each element contains the neighbors of the corresponding vertex
#'
#' @examples
#' \dontrun{
#' # Create a neighbor graph
#' g <- make_ring(5)
#' ng <- neighbor_graph(g)
#' 
#' # Get neighbors of vertex 1
#' n1 <- neighbors(ng, 1)
#' 
#' # Get all neighbors
#' all_n <- neighbors(ng)
#' }
#'
#' @seealso
#' * \code{\link{non_neighbors.neighbor_graph}} for finding non-neighboring vertices
#' * \code{\link{edges.neighbor_graph}} for getting all edges
#'
#' @export
neighbors.neighbor_graph <- function(x, i) {
  if (!missing(i)) {
    igraph::adjacent_vertices(x$G,i)
  } else {
    igraph::adjacent_vertices(x$G, 1:igraph::vcount(x$G))
  }

}

#' Node Density for neighbor_graph Objects
#'
#' @description
#' Computes the local density around each node based on the distances to its neighbors
#' in the feature space.
#'
#' @details
#' For each node i, computes the average squared Euclidean distance to its neighbors,
#' normalized by the square of the number of neighbors. This provides a measure of
#' how densely packed the neighborhood is in the feature space.
#'
#' The density for node i is computed as:
#' d_i = 1/(k_i^2) * sum((x_j - x_i)^2) for all neighbors j
#' where:
#' * k_i is the number of neighbors of node i
#' * x_i is the feature vector of node i
#' * x_j are the feature vectors of the neighbors
#'
#' @param x A neighbor_graph object
#' @param X A numeric matrix where rows represent data points and columns represent features
#'
#' @return A numeric vector containing the computed density for each node
#'
#' @examples
#' \dontrun{
#' # Create data matrix
#' X <- matrix(rnorm(50), 10, 5)
#' 
#' # Create neighbor graph
#' g <- make_ring(10)
#' ng <- neighbor_graph(g)
#' 
#' # Compute node densities
#' densities <- node_density(ng, X)
#' }
#'
#' @seealso
#' * \code{\link{neighbors.neighbor_graph}} for finding neighboring vertices
#'
#' @export
node_density.neighbor_graph <- function(x, X) {
  L <- neighbors(x)
  d <- sapply(seq_along(L), function(i) {
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
#' @description
#' Retrieves all edges from a neighbor_graph object as a two-column matrix.
#'
#' @details
#' Returns a matrix where each row represents an edge, with the first column
#' containing the source vertex index and the second column containing the
#' target vertex index.
#'
#' @param x A neighbor_graph object
#'
#' @return A numeric matrix with two columns:
#'   * Column 1: Source vertex indices
#'   * Column 2: Target vertex indices
#'
#' @examples
#' \dontrun{
#' # Create a neighbor graph
#' g <- make_ring(5)
#' ng <- neighbor_graph(g)
#' 
#' # Get all edges
#' edges <- edges(ng)
#' }
#'
#' @seealso
#' * \code{\link{neighbors.neighbor_graph}} for finding neighboring vertices
#'
#' @importFrom igraph as_data_frame
#' @export
edges.neighbor_graph <- function(x) {
  ret <- igraph::as_data_frame(x$G, what="edges")
  cbind(ret[,1], ret[,2])
}


#' Non-neighbors for neighbor_graph Objects
#'
#' @description
#' Identifies vertices that are not connected to a specified vertex in the graph.
#'
#' @details
#' For a given vertex i, returns all vertices that do not share an edge with i.
#' This is useful for finding disconnected components or analyzing graph structure.
#'
#' @param x A neighbor_graph object
#' @param i Integer specifying the vertex index for which to find non-neighbors
#'
#' @return A numeric vector containing indices of vertices that are not neighbors
#'   of the specified vertex i
#'
#' @examples
#' \dontrun{
#' # Create a neighbor graph
#' g <- make_ring(5)
#' ng <- neighbor_graph(g)
#' 
#' # Find non-neighbors of vertex 1
#' non_n <- non_neighbors(ng, 1)
#' }
#'
#' @seealso
#' * \code{\link{neighbors.neighbor_graph}} for finding neighboring vertices
#'
#' @importFrom igraph vcount
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

#' Print method for neighbor_graph objects
#'
#' @param x A neighbor_graph object
#' @param ... Additional arguments (not used)
#' @importFrom crayon bold blue green yellow white red
#' @export
print.neighbor_graph <- function(x, ...) {
  # Header
  cat("\n", crayon::bold(crayon::blue("Neighbor Graph Summary")), "\n")
  cat(crayon::white("─" %>% strrep(50)), "\n\n")
  
  # Basic information
  adj <- adjacency(x)
  n_vertices <- nrow(adj)
  n_edges <- sum(adj != 0) / 2  # Divide by 2 for undirected graphs
  sparsity <- 100 * (1 - n_edges / (n_vertices * (n_vertices - 1) / 2))
  
  cat(crayon::bold("📊 Basic Information:"), "\n")
  cat("   • Vertices:", crayon::green(n_vertices), "\n")
  cat("   • Edges:", crayon::green(sprintf("%.0f", n_edges)), "\n")
  cat("   • Sparsity:", crayon::green(sprintf("%.1f%%", sparsity)), "\n")
  cat("\n")
  
  # Weight information if weighted
  if ("weight" %in% names(x$params)) {
    weights <- adj@x[adj@x != 0]
    cat(crayon::bold("⚖️  Weight Distribution:"), "\n")
    cat("   • Range:", crayon::green(sprintf("[%.2f, %.2f]", min(weights), max(weights))), "\n")
    cat("   • Mean:", crayon::green(sprintf("%.2f", mean(weights))), "\n")
    cat("   • Median:", crayon::green(sprintf("%.2f", median(weights))), "\n")
    cat("\n")
  }
  
  # Graph parameters
  cat(crayon::bold("🔧 Parameters:"), "\n")
  for (param_name in names(x$params)) {
    param_value <- x$params[[param_name]]
    if (is.numeric(param_value)) {
      cat("   •", crayon::blue(param_name), ":", crayon::green(sprintf("%.2f", param_value)), "\n")
    } else {
      cat("   •", crayon::blue(param_name), ":", crayon::green(param_value), "\n")
    }
  }
  cat("\n")
  
  # Connectivity information
  degrees <- rowSums(adj != 0)
  avg_degree <- mean(degrees)
  max_degree <- max(degrees)
  isolated <- sum(degrees == 0)
  
  cat(crayon::bold("🔗 Connectivity:"), "\n")
  cat("   • Average degree:", crayon::green(sprintf("%.1f", avg_degree)), "\n")
  cat("   • Maximum degree:", crayon::green(max_degree), "\n")
  if (isolated > 0) {
    cat("   • Isolated vertices:", crayon::red(isolated), "\n")
  }
  cat("\n")
  
  # Footer
  cat(crayon::white("─" %>% strrep(50)), "\n")
  cat(crayon::bold(crayon::blue("Use str() for more details")), "\n\n")
  
  invisible(x)
}
