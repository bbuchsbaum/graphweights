

#' compute graph laplacian of a weight matrix
#'
#'
#' @param x the weight matrix
#' @param ... extra args
laplacian <- function(x,...) UseMethod("laplacian")

nvertices <- function(x, ...) UseMethod("nvertices")

#' neighbors of a set of nodes
neighbors <- function(x, i, ...) UseMethod("neighbors")

#' get weights of all neighbors of a set of nodes
adjacency <- function(x,...) UseMethod("adjacency")

#' get indices of the non-neighbors of a set of nodes
non_neighbors<- function(x,...) UseMethod("non_neighbors")

node_density <- function(x, ...) UseMethod("node_density")

class_means <- function(x, i, ...) UseMethod("class_means")

nclasses <- function(x) UseMethod("nclasses")

interclass_density <- function(x, X,...) UseMethod("interclass_density")

intraclass_density <- function(x, X,...) UseMethod("intraclass_density")

within_class_neighbors <- function(x, X, ...) UseMethod("within_class_neighbors")

between_class_neighbors <- function(x, X,...) UseMethod("between_class_neighbors")

intraclass_pairs <- function(x,...) UseMethod("intraclass_pairs")

interclass_pairs <- function(x,...) UseMethod("interclass_pairs")
