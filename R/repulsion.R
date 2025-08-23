#' @keywords internal
repulse_weight <- function(x1,x2, sigma=10) {
  #num <- sum((x1 - x2)^2)
  num <- sum((x1 - x2)^2)
  denom <- sqrt(sum(x1^2)) + sqrt(sum(x2^2))

  1/(sigma + num/denom)
}

#' Create a Repulsion Graph
#'
#' @description
#' Constructs a "repulsion graph" derived from an input graph `W` and a class structure graph `cg`.
#' The resulting graph retains only the edges from `W` that connect nodes belonging to *different* classes
#' according to `cg`. Edges connecting nodes within the same class are removed (or "repulsed").
#'
#' @details
#' This function takes an existing graph `W` (representing similarities, connections, etc.)
#' and filters it based on class labels. The `class_graph` object `cg` provides a binary
#' adjacency matrix where `1` indicates nodes belong to the *same* class. By taking the
#' complement (`!adjacency(cg)`), we get a mask where `1` indicates nodes belong to *different*
#' classes. Element-wise multiplication (`W * !adjacency(cg)`) effectively removes
#' within-class edges from `W`.
#'
#' The `method` argument controls whether the remaining edge weights are kept as is (`"weighted"`)
#' or converted to binary indicators (`"binary"`).
#'
#' This type of graph is useful when focusing on interactions *between* distinct groups or
#' classes within a larger network or dataset.
#'
#' @param W An input graph object. Can be a `Matrix` object (e.g., `dgCMatrix`) representing
#'   the adjacency matrix, or a `neighbor_graph` object. Edge weights are used if `method="weighted"`.
#' @param cg A `class_graph` object defining the class membership of the nodes in `W`.
#'   Must have the same dimensions as `W`.
#' @param method `character`. Specifies how to handle the weights of the remaining (between-class) edges:
#'   - `"weighted"` (default): Retains the original weights from `W`.
#'   - `"binary"`: Sets the weight of all remaining edges to 1.
#' @param threshold `numeric`. A threshold applied to the input graph `W` *before* filtering by class.
#'   Edges in `W` with weights strictly below this value are discarded. Default is 0.
#' @param norm_fac `numeric`. A normalization factor applied *only* if `method = "weighted"`.
#'   The weights of the retained edges are divided by this factor. Default is 1 (no normalization).
#'
#' @return A `repulsion_graph` object (inheriting from `neighbor_graph`), representing the filtered graph
#'   containing only between-class edges.
#'
#' @examples
#' library(Matrix)
#'
#' # --- Example 1: Simple Random Data ---
#' set.seed(123)
#' N <- 50
#' # Create a base graph (e.g., from features)
#' X <- matrix(rnorm(N * 5), N, 5)
#' W_adj <- Matrix(rsparsematrix(N, N, 0.1, symmetric = TRUE)) # Base adjacency
#' diag(W_adj) <- 0
#' W_ng <- neighbor_graph(W_adj) # Convert to neighbor_graph
#'
#' # Define class labels
#' labels <- factor(sample(1:3, N, replace = TRUE))
#' cg <- class_graph(labels)
#'
#' # Create a weighted repulsion graph
#' R_weighted <- repulsion_graph(W_ng, cg, method = "weighted")
#' print(R_weighted)
#' plot(R_weighted$G, vertex.color = labels, vertex.size=8, vertex.label=NA)
#' title("Weighted Repulsion Graph (Edges only between classes)")
#'
#' # Create a binary repulsion graph
#' R_binary <- repulsion_graph(W_ng, cg, method = "binary")
#' print(R_binary)
#'
#' # --- Example 2: Iris Dataset ---
#' data(iris)
#' X_iris <- as.matrix(iris[, 1:4])
#' labels_iris <- iris[, 5]
#' cg_iris <- class_graph(labels_iris)
#'
#' # Build a kNN graph first
#' W_iris_knn <- graph_weights(X_iris, k = 5, weight_mode = "heat", sigma = 0.7)
#'
#' # Create repulsion graph (weighted)
#' R_iris <- repulsion_graph(W_iris_knn, cg_iris, method = "weighted")
#' print(R_iris)
#' # Edges should primarily connect different iris species
#' plot(R_iris$G, vertex.color = as.numeric(labels_iris), vertex.size=5, vertex.label=NA)
#' title("Iris Repulsion Graph (k=5, heat weights)")
#'
#' @export
#' @importFrom assertthat assert_that
repulsion_graph <- function(W, cg, method = c("weighted", "binary"), threshold = 0,
                            norm_fac = 1) {
  method <- match.arg(method)
  assert_that(inherits(cg, "class_graph"), msg = "`cg` must be a class_graph object.")

  # Extract adjacency matrix and apply threshold
  if (inherits(W, "neighbor_graph")) {
    W_adj <- adjacency(W)
  } else if (inherits(W, "Matrix") || is.matrix(W)) {
    # Ensure W is a sparse Matrix for efficiency
    W_adj <- as(W, "sparseMatrix")
  } else {
    stop("`W` must be a `Matrix`, `matrix`, or `neighbor_graph` object")
  }

  N <- nrow(W_adj)
  cg_adj <- adjacency(cg)
  if (nrow(cg_adj) != N || ncol(cg_adj) != N) {
      stop("Dimensions of W and class_graph cg must match.")
  }

  # Apply threshold to input graph W
  if (threshold > 0) {
      W_adj@x[W_adj@x < threshold] <- 0
      W_adj <- drop0(W_adj) # Remove explicit zeros
  }

  # Create the repulsion mask (1 for between-class, 0 for within-class)
  repulsion_mask <- !adjacency(cg)

  # Apply the mask
  R_adj <- W_adj * repulsion_mask
  R_adj <- drop0(R_adj) # Clean up zeros created by multiplication

  # Handle weighted vs binary output
  if (method == "binary") {
    # Set non-zero weights to 1
    R_adj@x[R_adj@x != 0] <- 1
  } else {
    # Apply normalization factor if weighted
    if (norm_fac != 1) {
      if (norm_fac == 0) stop("norm_fac cannot be zero.")
      R_adj <- R_adj / norm_fac
    }
  }

  # Store parameters used
  params <- list(
      method = method,
      threshold = threshold,
      norm_fac = if(method == "weighted") norm_fac else NA_real_,
      original_W_type = class(W)[1],
      original_cg_levels = levels(cg$labels) # Or however levels are stored
  )

  # Create and return the new repulsion_graph object
  new_repulsion_graph(R_adj, params = params)
}

#' Constructor for Repulsion Graph Objects (Internal)
#'
#' @param adjacency_matrix The final sparse adjacency matrix of the repulsion graph.
#' @param params A list of parameters used in its construction.
#' @param ... Additional attributes to store (inherited from neighbor_graph).
#' @return An object of class `c("repulsion_graph", "neighbor_graph")`.
#' @keywords internal
new_repulsion_graph <- function(adjacency_matrix, params = list(), ...) {
    # Use the neighbor_graph constructor for Matrix input
    ng <- neighbor_graph.Matrix(adjacency_matrix, params = params, ...)
    # Add the specific class
    class(ng) <- c("repulsion_graph", class(ng))
    ng
}

#' Print method for repulsion_graph objects
#'
#' @param x A repulsion_graph object
#' @param ... Additional arguments passed to print
#' @importFrom crayon bold green blue red magenta
#' @export
print.repulsion_graph <- function(x, ...) {
    cat(bold(magenta("Repulsion Graph Object\n")))
    cat("----------------------\n")
    NextMethod("print") # Call the parent print method first

    # Add specific repulsion graph info
    cat(bold("Repulsion Params:\n"))
    cat("  Method:", blue(x$params$method), "\n")
    if (x$params$method == "weighted") {
        cat("  Normalization Factor:", blue(x$params$norm_fac %||% "NA"), "\n")
    }
    cat("  Input Threshold:", blue(x$params$threshold), "\n")
    cat("----------------------\n")
    invisible(x)
}

# Helper function for NULL default
`%||%` <- function(a, b) {
  if (is.null(a) || is.na(a)) b else a
}
