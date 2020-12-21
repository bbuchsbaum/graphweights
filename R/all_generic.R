

#' compute graph laplacian of a weight matrix
#'
#'
#' @param x the weight matrix
#' @param ... extra args
laplacian <- function(x,...) UseMethod("laplacian")


#' vertex count
#'
#' the number of vertices in the neighbor graph
#' @export
nvertices <- function(x, ...) UseMethod("nvertices")

#' neighbors of a set of nodes
#'
#' get the indices of neighbors of one or more vertices
#' @export
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

node_density <- function(x, i, ...) UseMethod("node_density")

class_means <- function(x, i, ...) UseMethod("class_means")

nclasses <- function(x) UseMethod("nclasses")

interclass_density <- function(x, X,...) UseMethod("interclass_density")

intraclass_density <- function(x, X,...) UseMethod("intraclass_density")

within_class_neighbors <- function(x, ng, ...) UseMethod("within_class_neighbors")

between_class_neighbors <- function(x, ng,...) UseMethod("between_class_neighbors")

intraclass_pairs <- function(x,...) UseMethod("intraclass_pairs")

interclass_pairs <- function(x,...) UseMethod("interclass_pairs")


edges <- function(x, ...) UseMethod("edges")

#' @export
search_result <- function(x, result, ...) UseMethod("search_result")

#' @export
dist_to_sim <- function(x, ...) UseMethod("dist_to_sim")

#' @export
find_nn <- function(x, ...) UseMethod("find_nn")

#' @export
find_nn_among <- function(x, ...) UseMethod("find_nn_among")

#' @export
find_nn_between <- function(x, ...) UseMethod("find_nn_between")

#' @export
neighbor_graph <- function(x, ...) UseMethod("neighbor_graph")


