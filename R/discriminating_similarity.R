#' Compute Discriminating Similarity
#'
#' Computes a similarity matrix that emphasizes discriminating features between classes.
#'
#' @param X Data matrix
#' @param labels Factor of class labels
#' @param k Number of nearest neighbors
#' @param sigma Bandwidth parameter for heat kernel
#' @return A similarity matrix emphasizing discriminating features
#' @export
discriminating_simililarity2 <- function(X, labels, k=nrow(X)/2, sigma=1) {
  # Compute initial similarity
  W <- graph_weights(X, k=k, weight_mode="euclidean", neighbor_mode="knn")
  
  # Get within-class and between-class matrices
  Ww <- create_label_matrix(labels, type="dense")
  Wb <- create_label_matrix(labels, type="dense", mode="distance")
  
  # Compute discriminating similarity
  D <- W * Ww / (W * Wb + .Machine$double.eps)
  
  # Apply heat kernel
  D <- heat_kernel(D, sigma=sigma)
  
  # Make symmetric
  (D + t(D))/2
}
