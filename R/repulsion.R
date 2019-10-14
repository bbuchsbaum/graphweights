
#' @param labels a set of labels defining the classes
#' @param W the neighborhood graph (weighted or binary)
#' @export
repulsion_graph <- function(labels, W) {
  cg <- class_graph(labels)
  acg <- 1-cg
  acg * W
}
