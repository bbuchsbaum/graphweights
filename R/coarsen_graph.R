#' Heavy Edge Potential Function
#'
#' This function assigns a potential to each edge based on its weight.
#'
#' @param edges data.table, the list of edges.
#' @return numeric vector, the potential of each edge.
#' @export
heavy_edge_potential <- function(edges) {
  return(edges$p_ij)  # Assume the potential is stored in the p_ij column
}

#' Normalize the Coarsening Matrix
#'
#' This function normalizes the coarsening matrix.
#'
#' @param C sparseMatrix, the coarsening matrix to be normalized.
#' @return sparseMatrix, the normalized coarsening matrix.
#' @import Matrix
#' @noRd
normalize_coarsening_matrix <- function(C) {
  norm_factors <- sqrt(rowSums(C^2))
  C <- Diagonal(x = 1 / norm_factors) %*% C
  return(C)
}

#' Randomized Edge Contraction Algorithm
#'
#' @param G sparseMatrix, the original graph as a sparse adjacency matrix.
#' @param T integer, the number of iterations.
#' @param potential_func function, a function that assigns a potential to each edge.
#' @return list, containing the coarsened graph, the contraction matrix, and the node mapping.
#' @import Matrix
#' @importFrom data.table data.table
#' @export
#' @examples
#' library(Matrix)
#' G <- sparseMatrix(i = c(1, 2, 3, 4), j = c(2, 3, 4, 1), x = c(1, 1, 1, 1))
#' result <- randomized_edge_contraction(G, 10, heavy_edge_potential)
#' coarsened_graph <- result$Gc
#' contraction_matrix <- result$C
#' node_mapping <- result$mapping
randomized_edge_contraction <- function(G, T, potential_func) {
  Gc <- G  # Initial contraction set is the graph itself
  num_vertices <- nrow(G)
  contraction_matrix <- Diagonal(num_vertices)
  active_nodes <- rep(TRUE, num_vertices)
  mapping <- seq_len(num_vertices)  # Initialize mapping of original nodes to coarsened nodes
  edges <- get_edges(Gc)

  for (t in 1:T) {
    if (nrow(edges) == 0) break  # Stop if no edges left

    phi <- edges$p_ij
    chosen_index <- sample(nrow(edges), 1, prob = phi)
    selected_edge <- edges[chosen_index, ]
    i <- selected_edge$i
    j <- selected_edge$j

    edges <- remove_connected_edges(edges, i, j)

    # Contract the edge
    Gc <- contract_edge(Gc, i, j)
    contraction_matrix <- update_coarsening_matrix(contraction_matrix, i, j)
    active_nodes[j] <- FALSE
    mapping[mapping == j] <- i  # Update mapping
  }

  # Create the coarsened graph with active nodes
  coarsened_nodes <- unique(mapping[active_nodes])
  Gc_coarsened <- Gc[coarsened_nodes, coarsened_nodes]

  # Create the coarsening matrix with N columns and n rows
  C_coarsened <- Matrix::sparseMatrix(
    i = rep(seq_along(coarsened_nodes), each = num_vertices),
    j = rep(seq_len(num_vertices), times = length(coarsened_nodes)),
    x = as.vector(contraction_matrix[, active_nodes]),
    dims = c(length(coarsened_nodes), num_vertices)
  )

  C_coarsened <- normalize_coarsening_matrix(C_coarsened)

  return(list(Gc = Gc_coarsened, C = C_coarsened, mapping = mapping))
}

#' Contract an Edge Between Two Nodes
#'
#' Modifies the adjacency matrix by merging node j into node i. This effectively removes node j from the graph.
#'
#' @param G sparseMatrix, the graph.
#' @param i Index of the node into which another node is merged.
#' @param j Index of the node that is merged into node i and removed.
#' @return The updated adjacency matrix after contracting the edge.
#' @noRd
contract_edge <- function(G, i, j) {
  # Merge j into i
  G[i, ] <- G[i, ] + G[j, ]
  G[, i] <- G[, i] + G[, j]

  # Avoid self-loop at i if i was connected to j
  G[i, i] <- 0

  # Clear the connections from and to j
  G[j, ] <- 0
  G[, j] <- 0

  return(G)
}

#' Update the Coarsening Matrix
#'
#' This function updates the coarsening matrix to reflect the contraction of two nodes.
#'
#' @param C sparseMatrix, the current coarsening matrix.
#' @param i integer, the node into which another node is merged.
#' @param j integer, the node that is merged into node i.
#' @return sparseMatrix, the updated coarsening matrix.
#' @noRd
update_coarsening_matrix <- function(C, i, j) {
  if (i != j) {
    C[, i] <- C[, i] + C[, j]
    C[, j] <- 0  # Zero out the column for j
  }
  return(C)
}

#' Heavy Edge Potential Function
#'
#' This function assigns a potential to each edge based on its weight.
#'
#' @param edges data.table, the list of edges.
#' @return numeric vector, the potential of each edge.
#' @export
heavy_edge_potential <- function(edges) {
  return(edges$p_ij)  # Assume the potential is stored in the p_ij column
}

#' Normalize the Coarsening Matrix
#'
#' This function normalizes the coarsening matrix.
#'
#' @param C sparseMatrix, the coarsening matrix to be normalized.
#' @return sparseMatrix, the normalized coarsening matrix.
#' @import Matrix
#' @noRd
normalize_coarsening_matrix <- function(C) {
  norm_factors <- sqrt(rowSums(C^2))
  C <- Diagonal(x = 1 / norm_factors) %*% C
  return(C)
}


#' @noRd
normalize_coarsening_matrix <- function(C) {
  # Normalize columns to sum to 1 (or other criteria depending on the application)
  col_sums <- sqrt(colSums(C^2))
  C <- C / col_sums
  return(C)
}




library(Matrix)

# Calculate the combinatorial Laplacian from an adjacency matrix
calculateLaplacian <- function(A) {
  if (!is(A, "Matrix")) stop("Input must be a sparse matrix.")
  if (nrow(A) != ncol(A)) stop("Adjacency matrix must be square.")

  D <- Diagonal(x = rowSums(A))
  L <- D - A
  return(L)
}

# Compute a matching based on the heavy edge rule
computeHeavyEdgeMatching <- function(A) {
  degree <- rowSums(A)
  edges <- which(A != 0, arr.ind = TRUE)
  weights <- A[edges] / pmax(degree[edges[, 1]], degree[edges[, 2]])

  matching <- data.frame(edges, weight = weights)
  matching <- matching[order(matching$weight, decreasing = TRUE), ]
  return(matching)
}

# Function to contract vertices and update the projection matrix
contractVertices <- function(A, u, v, P) {
  # Combine the rows and columns of vertices u and v
  A[u, ] <- A[u, ] + A[v, ]
  A[, u] <- A[, u] + A[, v]

  # Set the contracted vertex rows and columns to zero
  A[v, ] <- 0
  A[, v] <- 0

  # Update the diagonal to ensure no self-loops
  A[u, u] <- 0

  # Remove the contracted vertex from the matrix
  A <- A[-v, -v]

  # Update the projection matrix
  # Add the contribution of vertex v's projection to vertex u
  P[u, ] <- P[u, ] + P[v, ]
  P <- P[-v, ]

  return(list(A = A, P = P))
}

# Single level coarsening process using heavy edge matching (Algorithm 2)
singleLevelCoarsening <- function(L, A, targetSize, sigmaPrime) {
  N <- nrow(A)
  marked <- rep(FALSE, N)
  sigmaSquared <- 0
  P <- Diagonal(n = N)  # Initialize the projection matrix as an identity matrix

  matching <- computeHeavyEdgeMatching(A)

  while (nrow(matching) > 0 && N > targetSize && sqrt(sigmaSquared) <= sigmaPrime) {
    C <- matching[1, ]
    if (!marked[C$row] && !marked[C$col]) {
      marked[C$row] <- TRUE
      marked[C$col] <- TRUE
      N <- N - 1
      sigmaSquared <- sigmaSquared + C$weight
      result <- contractVertices(A, C$row, C$col, P)
      A <- result$A
      P <- result$P
    }
    matching <- matching[-1, ]
  }

  L_ell <- t(P) %*% L %*% P
  sigma_ell <- sqrt(sigmaSquared)

  return(list(L_ell = L_ell, sigma_ell = sigma_ell, P = P))
}

# Multi-level coarsening function
multiLevelCoarsening <- function(A, targetSize, epsilonPrime) {
  if (targetSize >= nrow(A) || targetSize <= 0) stop("Invalid target size.")

  L_ell <- calculateLaplacian(A)
  epsilon_ell <- 0
  P <- Diagonal(nrow(A))  # Initialize the projection matrix as an identity matrix

  while (nrow(A) > targetSize && epsilon_ell < epsilonPrime) {
    sigmaPrime <- (1 + epsilonPrime) / (1 + epsilon_ell) - 1
    result <- singleLevelCoarsening(L_ell, A, targetSize, sigmaPrime)
    L_ell <- result$L_ell
    sigma_ell <- result$sigma_ell
    P <- result$P %*% P  # Update the overall projection matrix
    epsilon_ell <- (1 + epsilon_ell) * (1 + sigma_ell) - 1
  }

  return(list(L_coarsened = L_ell, P = P))
}

# # Example usage
# adjMatrix <- Matrix(c(0, 1, 1, 0, 1, 0, 1, 1, 0), nrow = 3, sparse = TRUE)
# targetSize := 2
# result <- multiLevelCoarsening(adjMatrix, targetSize, epsilonPrime = 0.1)
# print(result$L_coarsened)
# print(result$P)
