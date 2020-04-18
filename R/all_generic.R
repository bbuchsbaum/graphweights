nrow          <- function(x) UseMethod("nrow")
nrow.default  <- base::nrow

ncol          <- function(x) UseMethod("ncol")
ncol.default  <- base::nrow


#' compute graph laplacian of a weight matrix
#'
#'
#' @param x the weight matrix
#' @param ... extra args
laplacian <- function(x,...) UseMethod("laplacian")


#' neighbors of a set of nodes
neighbors <- function(x, i, ...) UseMethod("neighbors")

#' get weights of all neighbors of a set of nodes
neighbor_weights <- function(x,...) UseMethod("neighbor_weights")

#' get indices of the non-neighbors of a set of nodes
non_neighbors<- function(x,...) UseMethod("non_neighbors")

class_means <- function(x, i, ...) UseMethod("class_means")

nclasses <- function(x) UseMethod("nclasses")

interclass_density <- function(x, X) UseMethod("interclass_density")

intraclass_density <- function(x, X) UseMethod("intraclass_density")

within_class_neighbors <- function(x, X, ...) UseMethod("within_class_neighbors")

between_class_neighbors <- function(x, X,...) UseMethod("between_class_neighbors")

