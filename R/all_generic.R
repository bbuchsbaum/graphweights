

#' Compute Graph Laplacian of a Weight Matrix
#'
#' This function computes the graph Laplacian of a given weight matrix. The graph Laplacian is defined as the difference between the degree matrix and the adjacency matrix.
#'
#' @param x The weight matrix representing the graph structure.
#' @param ... Additional arguments to be passed to specific implementations of the laplacian method.
#'
#' @return The graph Laplacian matrix of the given weight matrix.
#'
#' @export
laplacian <- function(x,...) UseMethod("laplacian")


#' Number of Vertices in Graph-like Objects
#'
#' Retrieve the number of vertices in a neighbor_graph object.
#'
#' @param x an object with a neighborhood
#' @param ... Additional arguments (currently ignored).
#'
#' @return The number of vertices in the neighbor_graph object.
#'
#' @examples
#' # Create a small example graph
#' adj_matrix <- matrix(c(0, 1, 1, 0, 1, 0, 1, 0, 0), nrow = 3, byrow = TRUE)
#' ng <- neighbor_graph(adj_matrix)
#' nvertices(ng) # Should return 3
#'
#' @export
nvertices <- function(x, ...) UseMethod("nvertices")

#' Neighbors of a Set of Nodes
#'
#' This function retrieves the indices of neighbors of one or more vertices in a given graph or graph-like object.
#'
#' @param x The graph or graph-like object in which to find the neighbors.
#' @param i The vertex or vertices for which to find the neighbors. Can be a single vertex index or a vector of vertex indices.
#' @param ... Additional arguments to be passed to specific implementations of the neighbors method.
#'
#' @return A list of vertex indices representing the neighbors of the specified vertices. The length of the list is equal to the number of input vertices, and each element in the list contains the neighbor indices for the corresponding input vertex.
#'
#' @examples
#' # Create an example graph
#' library(igraph)
#' g <- neighbor_graph(make_ring(5))
#'
#' # Get neighbors of vertex 1
#' n <- neighbors(g, 1)
neighbors <- function(x, i, ...) UseMethod("neighbors")

#' get weights of all neighbors of a set of nodes
#'
#' @param x the neighbor graph
#' @export
adjacency <- function(x,...) UseMethod("adjacency")


#' get indices of the non-neighbors of a set of nodes
#'
#' @param x the neighbor graph
#' @param ... extra args
#' @export
non_neighbors<- function(x,...) UseMethod("non_neighbors")

# node_density <- function(x, i, ...) UseMethod("node_density")

#' Class Means
#'
#' A generic function to compute the mean of each class.
#'
#' @param x An object.
#' @param i The index of the class for which the mean is to be computed.
#' @param ... Additional arguments passed to specific methods.
#'
#' @return A matrix or data frame representing the means of each class, the structure of which depends on the input object's class.
#'
#' @export
class_means <- function(x, i, ...) UseMethod("class_means")

#' Number of Classes
#'
#' A generic function to compute the number of classes in a graph.
#'
#' @param x An object.
#'
#' @return The number of classes in the input object.
#'
#' @export
nclasses <- function(x) UseMethod("nclasses")

# interclass_density <- function(x, X,...) UseMethod("interclass_density")
#
# intraclass_density <- function(x, X,...) UseMethod("intraclass_density")

#' Within-Class Neighbors
#'
#' A generic function to compute the within-class neighbors of a graph.
#'
#' @param x An object.
#' @param ng A neighbor graph object.
#' @param ... Additional arguments passed to specific methods.
#'
#' @return An object representing the within-class neighbors of the input graph, the structure of which depends on the input object's class.
#'
#' @export
within_class_neighbors <- function(x, ng, ...) UseMethod("within_class_neighbors")


#' Between-Class Neighbors
#'
#' A generic function to compute the between-class neighbors of a graph.
#'
#' @param x An object.
#' @param ng A neighbor graph object.
#' @param ... Additional arguments passed to specific methods.
#'
#' @return An object representing the between-class neighbors of the input graph, the structure of which depends on the input object's class.
#'
#' @export
between_class_neighbors <- function(x, ng,...) UseMethod("between_class_neighbors")

# intraclass_pairs <- function(x,...) UseMethod("intraclass_pairs")
# interclass_pairs <- function(x,...) UseMethod("interclass_pairs")

#' Edges for Graph-Like Objects
#'
#' Retrieve the edges of a graph-like object.
#'
#' @param x A graph-like object.
#' @param ... Further arguments passed to or from other methods.
#'
#' @return A matrix containing the edges of the graph-like object.
#' @export
edges <- function(x, ...) UseMethod("edges")

#' Neighbor Graph
#'
#' A generic function to create a neighbor_graph object.
#'
#' @param x An object.
#' @param ... Additional arguments passed to specific methods.
#'
#' @return A neighbor_graph object, the structure of which depends on the input object's class.
#'
#' @export
neighbor_graph <- function(x, ...) UseMethod("neighbor_graph")

#' Search result for nearest neighbor search
#'
#' @param x An object of class "nnsearcher".
#' @param result The result from the nearest neighbor search.
#'
#' @return An object with the class "nn_search".
#' @export
search_result <- function(x, result) {
  UseMethod("search_result")
}

#' Convert distance to similarity
#'
#' @param x An object representing a distance matrix.
#' @param ... Further arguments passed to or from other methods.
#'
#' @return A similarity matrix.
#' @export
dist_to_sim <- function(x, ...) {
  UseMethod("dist_to_sim")
}

#' Create an adjacency matrix
#'
#' @param x An object representing nearest neighbors.
#' @param ... Further arguments passed to or from other methods.
#'
#' @return An adjacency matrix.
#' @export
adjacency <- function(x, ...) {
  UseMethod("adjacency")
}

#' Find nearest neighbors
#'
#' @param x An object of class "nnsearcher".
#' @param ... Further arguments passed to or from other methods.
#'
#' @return A nearest neighbors result object.
#' @export
find_nn <- function(x, ...) {
  UseMethod("find_nn")
}

#' Find nearest neighbors among a subset
#'
#' @param x An object of class "nnsearcher" or "class_graph".
#' @param ... Further arguments passed to or from other methods.
#'
#' @return A nearest neighbors result object.
#' @export
find_nn_among <- function(x, ...) {
  UseMethod("find_nn_among")
}

#' Find nearest neighbors between two sets of data points
#'
#' @param x An object of class "nnsearcher".
#' @param ... Further arguments passed to or from other methods.
#'
#' @return A nearest neighbors result object.
#' @export
find_nn_between <- function(x, ...) {
  UseMethod("find_nn_between")
}


