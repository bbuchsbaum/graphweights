# Simplified ZoomOut Diagnostic Tests
# Focus on cases that should work reliably

library(Matrix)

# Load the zoomout function
source("../R/zoomout.R")

# Test 1: Larger Path Graph
test_larger_path <- function() {
  cat("=== Test 1: Larger Path Graph ===\n")
  
  # Create path graph: 1-2-3-4-5-6-7-8
  n <- 8
  W <- Matrix(0, n, n, sparse = TRUE)
  for(i in 1:(n-1)) {
    W[i,i+1] <- W[i+1,i] <- 1
  }
  
  # Simple identity functional map
  k0 <- 2
  C0 <- diag(k0)
  
  cat("Path graph (8 nodes), Identity C0\n")
  
  result <- zoomout(W, W, C0, k_max = 4,
                   normalized_laplacian = FALSE,
                   ncv = 10, seed = 42, verbose = TRUE)
  
  cat("C_final:\n")
  print(round(result$C_final, 3))
  cat("T_final:", result$T_final, "\n")
  
  # Check if reasonable
  is_permutation <- all(result$T_final >= 1 & result$T_final <= n)
  cat("Valid node indices?", is_permutation, "\n\n")
  
  return(result)
}

# Test 2: Grid Graph 
test_grid <- function() {
  cat("=== Test 2: Grid Graph (3x3) ===\n")
  
  # 3x3 grid graph
  # 1-2-3
  # | | |
  # 4-5-6  
  # | | |
  # 7-8-9
  W <- Matrix(0, 9, 9, sparse = TRUE)
  
  # Horizontal connections
  W[1,2] <- W[2,1] <- 1; W[2,3] <- W[3,2] <- 1
  W[4,5] <- W[5,4] <- 1; W[5,6] <- W[6,5] <- 1  
  W[7,8] <- W[8,7] <- 1; W[8,9] <- W[9,8] <- 1
  
  # Vertical connections
  W[1,4] <- W[4,1] <- 1; W[2,5] <- W[5,2] <- 1; W[3,6] <- W[6,3] <- 1
  W[4,7] <- W[7,4] <- 1; W[5,8] <- W[8,5] <- 1; W[6,9] <- W[9,6] <- 1
  
  C0 <- diag(3)
  
  cat("3x3 Grid graph, Identity C0\n")
  
  result <- zoomout(W, W, C0, k_max = 5,
                   normalized_laplacian = FALSE,
                   ncv = 12, seed = 42, verbose = TRUE)
  
  cat("C_final:\n")
  print(round(result$C_final, 3))
  cat("T_final:", result$T_final, "\n")
  
  # Should be identity mapping
  identity_check <- all(result$T_final == 1:9)
  cat("Perfect identity mapping?", identity_check, "\n\n")
  
  return(result)
}

# Test 3: Very Simple Case - Just check if algorithm runs
test_minimal <- function() {
  cat("=== Test 3: Minimal Working Case ===\n")
  
  # 4-node path: 1-2-3-4
  W <- Matrix(0, 4, 4, sparse = TRUE)
  W[1,2] <- W[2,1] <- 1
  W[2,3] <- W[3,2] <- 1
  W[3,4] <- W[4,3] <- 1
  
  # 1x1 functional map (minimal)
  C0 <- matrix(1, 1, 1)
  
  cat("4-node path, 1x1 C0\n")
  
  result <- zoomout(W, W, C0, k_max = 2,
                   normalized_laplacian = FALSE,
                   ncv = 6, seed = 42, verbose = TRUE)
  
  cat("C_final:\n")
  print(result$C_final)
  cat("T_final:", result$T_final, "\n")
  
  # Basic sanity checks
  valid_map <- all(result$T_final >= 1 & result$T_final <= 4)
  cat("Valid mapping?", valid_map, "\n")
  
  no_zero_indices <- all(result$T_final != 0)
  cat("No zero indices?", no_zero_indices, "\n\n")
  
  return(result)
}

# Test 4: Different initial maps
test_different_init <- function() {
  cat("=== Test 4: Different Initial Maps ===\n")
  
  # Use same 4-node path
  W <- Matrix(0, 4, 4, sparse = TRUE)
  W[1,2] <- W[2,1] <- 1
  W[2,3] <- W[3,2] <- 1
  W[3,4] <- W[4,3] <- 1
  
  cat("Testing different initial maps on 4-node path:\n")
  
  # Test identity
  C0_id <- diag(2)
  cat("Identity C0:\n")
  result1 <- zoomout(W, W, C0_id, k_max = 2,
                    normalized_laplacian = FALSE,
                    ncv = 6, seed = 42, verbose = FALSE)
  cat("T_final:", result1$T_final, "\n")
  
  # Test permutation
  C0_perm <- matrix(c(0,1,1,0), 2, 2)
  cat("Permutation C0:\n")
  result2 <- zoomout(W, W, C0_perm, k_max = 2,
                    normalized_laplacian = FALSE,
                    ncv = 6, seed = 42, verbose = FALSE)
  cat("T_final:", result2$T_final, "\n")
  
  # Test random
  set.seed(123)
  C0_rand <- matrix(rnorm(4), 2, 2)
  cat("Random C0:\n")
  result3 <- zoomout(W, W, C0_rand, k_max = 2,
                    normalized_laplacian = FALSE,
                    ncv = 6, seed = 42, verbose = FALSE)
  cat("T_final:", result3$T_final, "\n\n")
}

# Run tests
run_simple_tests <- function() {
  cat("ZoomOut Simple Diagnostic Tests\n")
  cat("===============================\n\n")
  
  test_larger_path()
  test_grid() 
  test_minimal()
  test_different_init()
  
  cat("Simple diagnostic tests completed!\n")
}

# Execute
if(sys.nframe() == 0) {
  run_simple_tests()
}