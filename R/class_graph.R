#' Construct a Class Graph
#'
#' Creates a graph structure where edges connect members of the same class. This is useful
#' for various machine learning tasks where class relationships need to be preserved.
#'
#' @param labels A vector of class labels that will be converted to factors if not already.
#' @param sparse Logical; whether to use sparse matrices in the computation. Default is TRUE.
#'              Using sparse matrices can significantly reduce memory usage for large datasets.
#' @param self_connections Logical; whether to include self-connections (diagonal elements set to 1).
#'              Default is FALSE, following standard graph theory conventions.
#'
#' @return A `class_graph` object, which is a list containing:
#' \itemize{
#'   \item{adjacency}{A matrix representing the adjacency structure of the graph.}
#'   \item{params}{A list of parameters used in graph construction.}
#'   \item{labels}{A vector of class labels.}
#'   \item{class_indices}{A list of vectors containing indices for each class.}
#'   \item{class_freq}{A table of class frequencies.}
#'   \item{levels}{A vector of unique class labels.}
#'   \item{classes}{Character string "class_graph" indicating the object type.}
#' }
#'
#' @examples
#' # Create a simple class graph from iris dataset
#' data(iris)
#' labels <- iris$Species
#' 
#' # Without self-connections (default)
#' cg <- class_graph(labels)
#' 
#' # With self-connections
#' cg_self <- class_graph(labels, self_connections = TRUE)
#'
#' @importFrom Matrix sparseVector tcrossprod Matrix t
#' @export
class_graph <- function(labels, sparse=TRUE, self_connections=FALSE) {
  # Input validation
  if (is.null(labels) || length(labels) == 0) {
    stop("labels cannot be NULL or empty")
  }
  
  # Check for NA values
  if (any(is.na(labels))) {
    stop("labels cannot contain NA values")
  }
  
  labels <- as.factor(labels)
  n <- length(labels)
  
  # Create empty sparse adjacency matrix
  adj <- Matrix::Matrix(0, n, n, sparse=TRUE)
  
  # Fill in edges for each class
  for (lev in levels(labels)) {
    idx <- which(labels == lev)
    if (length(idx) > 1) {
      # Create complete subgraph for class members
      for (i in 1:(length(idx)-1)) {
        for (j in (i+1):length(idx)) {
          adj[idx[i], idx[j]] <- 1
          adj[idx[j], idx[i]] <- 1  # Make symmetric
        }
      }
    }
  }
  
  # Set diagonal according to self_connections parameter
  if (self_connections) {
    diag(adj) <- 1
  }
  
  # Create class_graph object
  ret <- list(
    adjacency = adj,
    params = list(
      weight_mode = "binary",
      neighbor_mode = "supervised",
      self_connections = self_connections
    ),
    labels = labels,
    class_indices = split(1:length(labels), labels),
    class_freq = table(labels),
    levels = levels(labels),
    classes = "class_graph"
  )
  
  # Set class
  class(ret) <- c("class_graph", "neighbor_graph")
  
  ret
}

#'@export 
adjacency.class_graph <- function(x) {
  x$adjacency
}


#' Get Number of Classes in a class_graph Object
#'
#' @param x A class_graph object.
#'
#' @return An integer indicating the number of unique classes in the graph.
#'
#' @examples
#' data(iris)
#' cg <- class_graph(iris$Species)
#' n <- nclasses(cg)  # Should return 3 for iris dataset
#'
#' @export
nclasses.class_graph <- function(x) {
  length(x$levels)
}


#' Compute Class Means for a class_graph Object
#'
#' @param x A class_graph object.
#' @param X A numeric matrix or data frame where rows represent observations and columns represent features.
#'
#' @return A matrix where each row contains the mean values for a class, with row names corresponding to class labels.
#'
#' @examples
#' # Using iris dataset
#' data(iris)
#' X <- as.matrix(iris[, 1:4])
#' cg <- class_graph(iris$Species)
#' means <- class_means(cg, X)
#'
#' @export
class_means.class_graph <- function(x, X) {
  ret <- do.call(rbind, lapply(x$class_indices, function(i) {
    colMeans(X[i,,drop=FALSE])
  }))

  row.names(ret) <- names(x$class_indices)
  ret
}


#' Within-Class Neighbors for class_graph Objects
#'
#' Computes a neighbor graph that only includes connections between vertices of the same class.
#' This is useful for analyzing class-specific network structures.
#'
#' @param x A class_graph object containing class membership information.
#' @param ng A neighbor_graph object containing the base network structure.
#' @param ... Additional arguments passed to methods (currently unused).
#'
#' @return A neighbor_graph object containing only within-class connections.
#'
#' @examples
#' # Create example data
#' data(iris)
#' X <- as.matrix(iris[, 1:4])
#' cg <- class_graph(iris$Species)
#' ng <- neighbor_graph(X)
#' within_ng <- within_class_neighbors(cg, ng)
#'
#' @export
within_class_neighbors.class_graph <- function(x, ng, ...) {
  Ac <- adjacency(x)
  An <- adjacency(ng)
  Aout <- Ac * An
  neighbor_graph(Aout)
}


#' Between-Class Neighbors for class_graph Objects
#'
#' Computes a neighbor graph that only includes connections between vertices of different classes.
#' This is useful for analyzing relationships between different classes in the network.
#'
#' @param x A class_graph object containing class membership information.
#' @param ng A neighbor_graph object containing the base network structure.
#' @param ... Additional arguments passed to methods (currently unused).
#'
#' @return A neighbor_graph object containing only between-class connections.
#'
#' @examples
#' # Create example data
#' data(iris)
#' X <- as.matrix(iris[, 1:4])
#' cg <- class_graph(iris$Species)
#' ng <- neighbor_graph(X)
#' between_ng <- between_class_neighbors(cg, ng)
#'
#' @export
between_class_neighbors.class_graph <- function(x, ng, ...) {
  Ac <- !adjacency(x)
  An <- adjacency(ng)
  Aout <- Ac * An
  neighbor_graph(Aout)
}


#' Compute Discriminating Distance for Similarity Graph
#'
#' Computes a discriminating distance matrix for a similarity graph based on class labels.
#' The function modifies the weights within and between classes to enhance class discrimination,
#' making it particularly useful for classification and clustering tasks.
#'
#' @param X A numeric matrix or data frame where rows are observations and columns are features.
#' @param k Integer; number of nearest neighbors to consider. Defaults to half the number of labels.
#' @param sigma Numeric; scaling factor for the heat kernel. If missing, it will be estimated from the data.
#' @param labels Factor or vector containing class labels for each observation.
#'
#' @return A numeric matrix representing the discriminating distances between observations.
#'
#' @examples
#' # Create example data
#' set.seed(123)
#' X <- matrix(rnorm(100 * 4), 100, 4)
#' labels <- factor(rep(1:5, each = 20))
#' 
#' # Compute discriminating distance
#' D <- discriminating_distance(X, k = 10, sigma = 0.7, labels = labels)
#'
#' @references
#' Cai, D., He, X., Zhou, K., Han, J., & Bao, H. (2007). 
#' Locality Sensitive Discriminant Analysis. 
#' In IJCAI (Vol. 2007, pp. 708-713).
#'
#' @export
discriminating_distance <- function(X, k=nrow(X)/2, sigma, labels) {
  # Compute initial similarity
  W <- graph_weights(X, k=k, weight_mode="euclidean", neighbor_mode="knn")
  
  # Get within-class and between-class matrices
  Ww <- create_label_matrix(labels, type="sparse")
  Wb <- create_label_matrix(labels, type="sparse", mode="distance")
  
  # Compute weighted similarities
  Ww2 <- W * Ww
  Wb2 <- W * Wb
  
  # Get indices
  wind <- which(Ww2 > 0)
  bind <- which(Wb2 > 0)
  
  # Apply heat kernel
  if (missing(sigma)) {
    sigma <- estimate_sigma(X, prop=.1)
  }

  hw <- inverse_heat_kernel(W[wind], sigma)
  hb <- inverse_heat_kernel(W[bind], sigma)
  
  # Create output matrix
  out <- Matrix::Matrix(0, nrow=nrow(W), ncol=ncol(W), sparse=TRUE)
  out[wind] <- hw * (1-hw)
  out[bind] <- hb * (1+hb)
  
  # Make symmetric
  (out + Matrix::t(out))/2
}

#' Compute Similarity Graph Weighted by Class Structure
#'
#' Computes a similarity graph that incorporates class structure information.
#' This method is particularly effective for pattern recognition tasks such as
#' face recognition and handwriting digit classification.
#'
#' @param X A numeric matrix or data frame where rows are observations and columns are features.
#' @param k Integer; number of nearest neighbors to consider. Defaults to half the number of labels.
#' @param sigma Numeric; scaling factor for the heat kernel.
#' @param cg A class_graph object containing class structure information.
#' @param threshold Numeric; threshold value for the class graph (default: 0.01).
#'
#' @return A numeric matrix representing the weighted similarity graph.
#'
#' @examples
#' # Create example data
#' set.seed(123)
#' X <- matrix(rnorm(100 * 4), 100, 4)
#' labels <- factor(rep(1:5, each = 20))
#' cg <- class_graph(labels)
#' 
#' # Compute discriminating similarity
#' W <- discriminating_similarity(X, k = 10, sigma = 0.7, cg = cg)
#'
#' @references
#' Zhang, W., Xue, X., Sun, Z., Guo, Y., & Lu, H. (2007).
#' Local similarity and diversity preserving discriminant projection for face
#' and handwriting digits recognition.
#' Neurocomputing, 70(7-9), 1454-1461.
#'
#' @export
discriminating_similarity <- function(X, k=length(labels)/2, sigma, cg, threshold=.01) {
  #Wknn <- graph_weights(X)

  Wall <- graph_weights(X, k=k, weight_mode="heat", neighbor_mode="knn",sigma=sigma)
  W_adj <- adjacency(Wall)  # Get adjacency matrix

  Ww <- create_label_matrix(labels, type="sparse")
  Wb <- create_label_matrix(labels, type="sparse", mode="distance")

  Ww2 <- W_adj * Ww
  Wb2 <- W_adj * Wb

  wind <- which(Ww2 > 0)
  bind <- which(Wb2 > 0)

  hw <- heat_kernel(W_adj[wind], sigma)
  hb <- heat_kernel(W_adj[bind], sigma)

  W_adj[wind] <- hw * (1+hw)
  W_adj[bind] <- hb * (1-hb)
  W_adj
}

#' Print method for class_graph objects
#'
#' @param x A class_graph object
#' @param ... Additional arguments passed to print
#'
#' @importFrom crayon bold blue green red yellow white
#' @export
print.class_graph <- function(x, ...) {
  # Header
  cat("\n", crayon::bold(crayon::blue("Class Graph Summary")), "\n")
  cat(crayon::white("─" %>% strrep(50)), "\n\n")
  
  # Basic information
  n_vertices <- if (!is.null(x$adjacency)) nrow(x$adjacency) else length(x$labels)
  
  cat(crayon::bold("📊 Basic Information:"), "\n")
  cat("   • Number of vertices:", crayon::green(n_vertices), "\n")
  cat("   • Number of classes:", crayon::green(length(x$levels)), "\n")
  
  # Only show sparsity if adjacency matrix exists
  if (!is.null(x$adjacency)) {
    sparsity <- 100 * sum(x$adjacency != 0) / (nrow(x$adjacency) * ncol(x$adjacency))
    cat("   • Sparsity:", crayon::green(sprintf("%.2f%%", sparsity)), "\n")
  }
  cat("\n")
  
  # Class distribution
  cat(crayon::bold("📈 Class Distribution:"), "\n")
  max_freq <- max(x$class_freq)
  bar_width <- 30
  
  for (i in seq_along(x$class_freq)) {
    class_name <- names(x$class_freq)[i]
    freq <- x$class_freq[i]
    bar_len <- round(freq / max_freq * bar_width)
    
    # Create bar with gradient colors
    bar <- crayon::yellow("█" %>% strrep(bar_len))
    
    # Format output
    cat(sprintf("   %-15s %4d %s\n", 
        crayon::blue(class_name), 
        freq,
        bar))
  }
  cat("\n")
  
  # Edge information
  cat(crayon::bold("🔗 Edge Information:"), "\n")
  
  if (is.null(x$adjacency)) {
    cat("   • No adjacency matrix found in graph\n")
  } else {
    # Count edges properly
    adj_mat <- as.matrix(x$adjacency)
    # Only count upper triangle since matrix is symmetric
    total_edges <- sum(adj_mat[upper.tri(adj_mat)] != 0)
    
    if (total_edges == 0) {
      cat("   • No edges found in the graph\n")
    } else {
      cat("   • Total edges:", crayon::green(sprintf("%d", total_edges)), "\n")
      
      # Within-class edges
      within_class_edges <- 0
      for (class_idx in x$class_indices) {
        if (length(class_idx) > 1) {  # Only process classes with more than one member
          class_mat <- adj_mat[class_idx, class_idx]
          within_class_edges <- within_class_edges + 
            sum(class_mat[upper.tri(class_mat)] != 0)
        }
      }
      
      between_class_edges <- total_edges - within_class_edges
      
      # Calculate percentages safely
      within_pct <- if (total_edges > 0) 100 * within_class_edges / total_edges else 0
      between_pct <- if (total_edges > 0) 100 * between_class_edges / total_edges else 0
      
      cat("   • Within-class edges:", crayon::green(sprintf("%d", within_class_edges)), 
          sprintf("(%.1f%%)", within_pct), "\n")
      cat("   • Between-class edges:", crayon::green(sprintf("%d", between_class_edges)),
          sprintf("(%.1f%%)", between_pct), "\n")
    }
  }
  cat("\n")
  
  # Footer
  cat(crayon::white("─" %>% strrep(50)), "\n")
  cat(crayon::bold(crayon::blue("Use str() for more details")), "\n\n")
  
  invisible(x)
}
