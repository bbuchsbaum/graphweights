# Diagnostic Tests for ZoomOut Implementation
# Simple cases where we can predict expected behavior

library(Matrix)
library(testthat)

# Load the zoomout function
source("../R/zoomout.R")

# Test 1: Identity Case - Same Graph, Identity Functional Map
test_identity_case <- function() {
  cat("=== Test 1: Identity Case ===\n")
  
  # Create simple path graph: 1-2-3-4-5
  n <- 5
  W <- Matrix(0, n, n, sparse = TRUE)
  W[1,2] <- W[2,1] <- 1
  W[2,3] <- W[3,2] <- 1  
  W[3,4] <- W[4,3] <- 1
  W[4,5] <- W[5,4] <- 1
  
  # Perfect identity functional map
  k0 <- 2
  C0 <- diag(k0)  # Identity matrix
  
  cat("Input: Path graph (5 nodes), Identity C0 =\n")
  print(C0)
  
  # Run ZoomOut (use small k_max due to graph size limits)
  result <- zoomout(W, W, C0, k_max = 2,  # Very small k_max for 5-node graph
                   normalized_laplacian = FALSE, 
                   ncv = 8,  # Explicit ncv setting
                   seed = 42, verbose = TRUE)
  
  cat("Final functional map C_final =\n")
  print(round(result$C_final, 4))
  
  cat("Final pointwise map T_final:", result$T_final, "\n")
  
  # Check if we get identity pointwise map
  expected_identity <- 1:n
  cat("Expected identity mapping:", expected_identity, "\n")
  cat("Matches identity?", all(result$T_final == expected_identity), "\n\n")
  
  return(result)
}

# Test 2: Simple Permutation - Reversed Graph
test_permutation_case <- function() {
  cat("=== Test 2: Simple Permutation Case ===\n")
  
  # Source: 1-2-3 (path)
  W_M <- Matrix(0, 3, 3, sparse = TRUE)
  W_M[1,2] <- W_M[2,1] <- 1
  W_M[2,3] <- W_M[3,2] <- 1
  
  # Target: 3-2-1 (reversed path - should be isomorphic)
  W_N <- W_M  # Same structure
  
  # Initial functional map suggesting reversal
  # For 2x2 case, try simple permutation matrix
  C0 <- matrix(c(0, 1, 1, 0), 2, 2)  # Simple permutation
  
  cat("Source graph: 1-2-3\n")
  cat("Target graph: 1-2-3 (same)\n") 
  cat("Initial C0 (permutation) =\n")
  print(C0)
  
  result <- zoomout(W_M, W_N, C0, k_max = 2,  # Small k_max
                   normalized_laplacian = FALSE,
                   seed = 42, verbose = TRUE)
  
  cat("Final C_final =\n")
  print(round(result$C_final, 4))
  cat("T_final:", result$T_final, "\n\n")
  
  return(result)
}

# Test 3: Degenerate Case - Single Node
test_single_node <- function() {
  cat("=== Test 3: Single Node Case ===\n")
  
  # Single isolated node
  W_single <- Matrix(0, 1, 1, sparse = TRUE)
  C0 <- matrix(1, 1, 1)  # 1x1 identity
  
  cat("Single node graph, C0 = 1x1 identity\n")
  
  tryCatch({
    result <- zoomout(W_single, W_single, C0, k_max = 1,
                     normalized_laplacian = FALSE, 
                     seed = 42, verbose = TRUE)
    cat("Result successful - T_final:", result$T_final, "\n")
  }, error = function(e) {
    cat("Expected error for single node:", e$message, "\n")
  })
  
  cat("\n")
}

# Test 4: Triangle Graph - Symmetric Case
test_triangle_case <- function() {
  cat("=== Test 4: Triangle Graph Case ===\n")
  
  # Complete triangle: all nodes connected
  W_tri <- Matrix(0, 3, 3, sparse = TRUE)
  W_tri[1,2] <- W_tri[2,1] <- 1
  W_tri[1,3] <- W_tri[3,1] <- 1  
  W_tri[2,3] <- W_tri[3,2] <- 1
  
  # Identity functional map
  C0 <- diag(2)
  
  cat("Complete triangle graph, Identity C0\n")
  
  result <- zoomout(W_tri, W_tri, C0, k_max = 2,  # Small k_max
                   normalized_laplacian = FALSE,
                   seed = 42, verbose = TRUE)
  
  cat("C_final =\n")
  print(round(result$C_final, 4))
  cat("T_final:", result$T_final, "\n")
  
  # For symmetric graph, should get some valid permutation
  is_permutation <- all(sort(result$T_final) == 1:3)
  cat("Is valid permutation?", is_permutation, "\n\n")
  
  return(result)
}

# Test 5: Star Graph - Hub Structure
test_star_case <- function() {
  cat("=== Test 5: Star Graph Case ===\n")
  
  # Star graph: node 1 connected to all others
  n <- 5
  W_star <- Matrix(0, n, n, sparse = TRUE)
  for(i in 2:n) {
    W_star[1,i] <- W_star[i,1] <- 1
  }
  
  C0 <- diag(2)
  
  cat("Star graph (5 nodes), center = node 1\n")
  
  result <- zoomout(W_star, W_star, C0, k_max = 2,  # Small k_max
                   normalized_laplacian = FALSE,
                   ncv = 8,
                   seed = 42, verbose = TRUE)
  
  cat("C_final =\n") 
  print(round(result$C_final, 4))
  cat("T_final:", result$T_final, "\n")
  
  # Node 1 should map to itself (center of star)
  center_maps_correctly <- result$T_final[1] == 1
  cat("Center node maps to itself?", center_maps_correctly, "\n\n")
  
  return(result)
}

# Test 6: Verify Mathematical Properties
test_mathematical_properties <- function(result) {
  cat("=== Test 6: Mathematical Properties ===\n")
  
  C <- result$C_final
  Phi_M <- result$Phi_M
  Phi_N <- result$Phi_N
  
  # Check if C has reasonable magnitude
  cat("||C||_F =", norm(C, "F"), "\n")
  
  # Check orthogonality properties if possible
  if(nrow(C) == ncol(C)) {
    orthogonality <- norm(t(C) %*% C - diag(nrow(C)), "F")
    cat("||C^T C - I||_F =", orthogonality, "\n")
  }
  
  # Check if eigenvectors are orthonormal
  if(ncol(Phi_M) > 1) {
    M_orthogonality <- norm(t(Phi_M) %*% Phi_M - diag(ncol(Phi_M)), "F")
    cat("||Phi_M^T Phi_M - I||_F =", M_orthogonality, "\n")
  }
  
  if(ncol(Phi_N) > 1) {
    N_orthogonality <- norm(t(Phi_N) %*% Phi_N - diag(ncol(Phi_N)), "F")  
    cat("||Phi_N^T Phi_N - I||_F =", N_orthogonality, "\n")
  }
  
  cat("\n")
}

# Run all tests
run_diagnostic_tests <- function() {
  cat("ZoomOut Diagnostic Tests\n")
  cat("========================\n\n")
  
  # Run basic tests
  result1 <- test_identity_case()
  result2 <- test_permutation_case() 
  test_single_node()
  result4 <- test_triangle_case()
  result5 <- test_star_case()
  
  # Test mathematical properties on one result
  test_mathematical_properties(result1)
  
  cat("All diagnostic tests completed!\n")
}

# Execute if run directly
if(sys.nframe() == 0) {
  run_diagnostic_tests()
}