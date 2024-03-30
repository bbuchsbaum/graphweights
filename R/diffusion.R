# Install and load the necessary packages

# Define the function to compute the Markov diffusion kernel
compute_markov_diffusion_kernel <- function(A, t) {
  # Convert the adjacency matrix A into a transition probability matrix P
  P <- A / rowSums(A)

  # Compute the matrix exponential for P^t
  P_t <- expm(P * t)

  # Compute Z(t) matrix
  Z <- (diag(nrow(P)) - P) ^ (-1) * (diag(nrow(P)) - P_t) * P / t

  # Compute the Markov diffusion kernel
  K_MD <- Z %*% t(Z)

  return(K_MD)
}

compute_markov_diffusion_kernel_eigen <- function(A, t) {
  # Convert the adjacency matrix A into a transition probability matrix P
  P <- A / rowSums(A)

  # Compute the eigenvectors (U) and eigenvalues (L) of the matrix P
  eigen_decomp <- eigen(P)
  U <- eigen_decomp$vectors
  L <- diag(eigen_decomp$values)

  # Compute the matrix power of L using the input parameter t
  L_t <- L ^ t

  # Compute the Markov diffusion kernel (K_MD) as a product of the eigenvectors and the matrix power of eigenvalues
  K_MD <- U %*% L_t %*% t(U)

  return(K_MD)
}

#
# # Create a large sparse adjacency matrix (replace this with your actual matrix)
# A <- rsparsematrix(1000, 1000, density = 0.01)
#
# # Set the number of steps for the diffusion process
# t <- 10
#
# # Compute the Markov diffusion kernel for the given adjacency matrix
# K_MD <- compute_markov_diffusion_kernel(A, t)
#
# # Show the result
# print(K_MD)
#
#
#
# # Create a large sparse adjacency matrix (replace this with your actual matrix)
# A <- rsparsematrix(1000, 1000, density = 0.01)
#
# # Set the number of steps for the diffusion process
# t <- 10
#
# # Compute the Markov diffusion kernel for the given adjacency matrix
# K_MD <- compute_markov_diffusion_kernel(A, t)
#
# # Show the result
# print(K_MD)
#
#
# library(Matrix)
# library(igraph)

# Function to compute the diffusion distance and diffusion map
compute_diffusion <- function(adj_matrix, t) {
  # Convert the adjacency matrix to a graph
  g <- graph_from_adjacency_matrix(adj_matrix, mode = "undirected", weighted = TRUE)

  # Compute the degree matrix
  deg_matrix <- diag(degree(g, mode = "all"))

  # Compute the Laplacian matrix
  laplacian_matrix <- deg_matrix - adj_matrix

  # Compute the normalized Laplacian matrix
  norm_laplacian_matrix <- diag(1/sqrt(diag(deg_matrix))) %*% laplacian_matrix %*% diag(1/sqrt(diag(deg_matrix)))

  # Compute the eigenvectors and eigenvalues of the normalized Laplacian matrix
  eigen_decomposition <- eigen(norm_laplacian_matrix)
  eigenvalues <- eigen_decomposition$values
  eigenvectors <- eigen_decomposition$vectors

  # Sort eigenvectors and eigenvalues by eigenvalues (ascending)
  sorted_indices <- order(eigenvalues)
  sorted_eigenvalues <- eigenvalues[sorted_indices]
  sorted_eigenvectors <- eigenvectors[, sorted_indices]

  # Compute the diffusion map
  diffusion_map <- t(sorted_eigenvectors) %*% diag(sorted_eigenvalues^t) %*% t(sorted_eigenvectors)

  # Compute the diffusion distance matrix
  diffusion_distance_matrix <- matrix(0, nrow = nrow(diffusion_map), ncol = nrow(diffusion_map))
  for (i in 1:nrow(diffusion_map)) {
    for (j in 1:nrow(diffusion_map)) {
      diffusion_distance_matrix[i, j] <- sum((diffusion_map[i, ] - diffusion_map[j, ])^2)
    }
  }

  list(diffusion_map = diffusion_map, diffusion_distance_matrix = diffusion_distance_matrix)
}

# # Example adjacency matrix
# adj_matrix <- matrix(c(0, 1, 0, 0,
#                        1, 0, 1, 0,
#                        0, 1, 0, 1,
#                        0, 0, 1, 0), nrow = 4, ncol = 4)
#
# # Set the time 't'
# t <- 1
#
# # Compute the diffusion distance and diffusion map
# result <- compute_diffusion(adj_matrix, t)
#
# # Print the diffusion map and diffusion distance matrix
# print(result$diffusion_map)
# print(result$diffusion_distance_matrix)
#
