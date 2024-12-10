library(testthat)
library(Matrix)

test_that("Coarsening a simple path graph preserves spectral properties", {
  # Create a simple path graph of 4 nodes:
  # 1 -- 2 -- 3 -- 4
  # Adjacency:
  # 1:2, 2:1,3; 3:2,4; 4:3
  W <- sparseMatrix(
    i = c(1,2,3),
    j = c(2,3,4),
    x = c(1,1,1),
    dims = c(4,4),
    symmetric = TRUE
  )
  
  # Coarsen the graph deterministically for testing
  result <- rec_coarsen(W, T=1, deterministic_first_edge=TRUE)

  # Check dimensions of the coarsened graph
  # After one contraction, we must have fewer or equal number of coarse vertices.
  # Given that one edge was contracted, we expect n < N.
  n <- nrow(result$W_c)
  expect_lt(n, 4)
  
  # Check that W_c is symmetric and non-negative
  expect_true(isSymmetric(result$W_c))
  expect_true(all(result$W_c@x >= 0))
  
  # Check coarsening matrix C dimensions: C is n x N
  expect_equal(dim(result$C), c(n, 4))
  
  # Check that each row of C has norm 1
  Cd <- as.matrix(result$C)
  row_norms <- rowSums(Cd^2)
  expect_true(all(abs(row_norms - 1) < 1e-12))
  
  # Since we know coarsening doesn't produce a trivial adjacency,
  # let's check if the coarsened Laplacian L_c is a good spectral approximation
  # for the first few eigenvalues.
  
  # Compute original Laplacian
  d <- rowSums(W)
  L <- Diagonal(x=d) - W
  # Compute first 2 eigenvalues of L
  orig_eigs <- sort(eigen(L, symmetric=TRUE, only.values=TRUE)$values)[1:2]
  
  # Compute first 2 eigenvalues of L_c
  L_c <- result$L_c
  coars_eigs <- sort(eigen(as.matrix(L_c), symmetric=TRUE, only.values=TRUE)$values)[1:2]
  
  # We do not expect identical eigenvalues, but a coarse approximation.
  # We'll use a very loose tolerance here, just checking they are of the same order.
  # This is a heuristic check. Depending on the graph and method, you may adjust the tolerance.
  expect_true(abs(orig_eigs[2] - coars_eigs[2]) < 0.8 * orig_eigs[2] + 1e-1)
})

test_that("Coarsening a star graph preserves dimension and basic structure", {
  # A star graph with 5 nodes: center=1 connected to (2,3,4,5)
  N <- 5
  W <- sparseMatrix(
    i = c(1,1,1,1),
    j = c(2,3,4,5),
    x = c(1,1,1,1),
    dims = c(N,N),
    symmetric = TRUE
  )

  # Coarsen once
  result <- rec_coarsen(W, T=1, deterministic_first_edge=TRUE)
  
  # After one contraction, we expect n < N
  n <- nrow(result$W_c)
  expect_lt(n, N)
  
  # Check symmetry and non-negativity
  expect_true(isSymmetric(result$W_c))
  expect_true(all(result$W_c@x >= 0))
  
  # Check dimensions of C
  expect_equal(dim(result$C)[1], n)
  expect_equal(dim(result$C)[2], N)
  
  # Check normalization of rows in C
  Cd <- as.matrix(result$C)
  expect_true(all(abs(rowSums(Cd^2) - 1) < 1e-12))
  
  # Again, spectral property check
  d <- rowSums(W)
  L <- Diagonal(x=d) - W
  orig_eigs <- sort(eigen(L, symmetric=TRUE, only.values=TRUE)$values)[1:2]
  
  L_c <- result$L_c
  coars_eigs <- sort(eigen(as.matrix(L_c), symmetric=TRUE, only.values=TRUE)$values)[1:2]
  
  # Loose check on spectral approximation
  expect_true(abs(orig_eigs[2] - coars_eigs[2]) < 0.8 * orig_eigs[2] + 1e-1)
})