library(testthat)
library(Matrix)

context("coarsen_graph tests")

test_that("Coarsening a simple path graph preserves basic properties", {
  # Create a simple path graph of 4 nodes:
  # 1 -- 2 -- 3 -- 4
  W <- sparseMatrix(
    i = c(1,2,3),
    j = c(2,3,4),
    x = c(1,1,1),
    dims = c(4,4),
    symmetric = TRUE
  )

  W <- as(W, "dgCMatrix")

  # Coarsen the graph
  # Choose k < N, say k=2 for a small subspace,
  # epsilon_prime=0.5 is arbitrary, n=2 means we want to reduce to 2 vertices or so.
  result <- coarsen_graph_multilevel(W, n=2, k=2, epsilon_prime=0.5, verbose=FALSE)
  
  A_c <- result$A_coarse
  P <- result$P
  info <- result$info
  
  # Check dimensions
  expect_lt(nrow(A_c), 4)  # we expect fewer vertices than original
  expect_lt(ncol(P), 5)    # P should map from N to n, so n < N
  
  # Check that A_c is symmetric and non-negative
  expect_true(isSymmetric(A_c))
  expect_true(all(A_c@x >= 0))
  
  # Check P dimensions: P is n' x N
  expect_equal(ncol(P), 4)
  expect_equal(nrow(P), nrow(A_c))
  
  # Check that rows of P have non-zero norm. For Laplacian consistent coarsening,
  # rows correspond to sets of contracted vertices.
  # Depending on the construction, the norm may not be exactly 1, but let's verify 
  # they are not degenerate and follow the rules (should not be zero).
  row_norms <- sqrt(rowSums(P^2))
  expect_true(all(row_norms > 1e-15))
  
  # Check spectral approximation (very loose):
  # Compute original Laplacian
  d <- rowSums(W)
  L <- Diagonal(x=d) - W
  orig_eigs <- sort(eigen(as.matrix(L), symmetric=TRUE, only.values=TRUE)$values)[1:2]
  
  # Compute coarse Laplacian
  d_c <- rowSums(A_c)
  L_c <- Diagonal(x=d_c) - A_c
  coars_eigs <- sort(eigen(as.matrix(L_c), symmetric=TRUE, only.values=TRUE)$values)[1:2]
  
  # Loose check: the second smallest eigenvalue (Fiedler value) should not diverge completely
  # We'll just ensure it's finite and not negative
  expect_true(is.finite(coars_eigs[2]))
  expect_true(coars_eigs[2] >= -1e-12)
  
  # A more direct spectral approximation check could involve restricted similarity,
  # but here we just ensure no pathological values.
})

test_that("Coarsening a star graph reduces dimensions and keeps structure", {
  # A star graph with 5 nodes: center=1 connected to (2,3,4,5)
  N <- 5
  W <- sparseMatrix(
    i = c(1,1,1,1),
    j = c(2,3,4,5),
    x = c(1,1,1,1),
    dims = c(N,N),
    symmetric = TRUE
  )
  W <- as(W, "dgCMatrix")
  
  # Coarsen with n=3 to reduce size, k=2 for subspace dimension
  result <- coarsen_graph_multilevel(W, n=3, k=2, epsilon_prime=0.5, verbose=FALSE)
  
  A_c <- result$A_coarse
  P <- result$P
  info <- result$info
  
  # After coarsening, expect fewer than N vertices
  expect_lt(nrow(A_c), N)
  
  # Check symmetry and non-negativity of A_c
  expect_true(isSymmetric(A_c))
  expect_true(all(A_c@x >= 0))
  
  # Check P dimensions: P is n' x N
  expect_equal(ncol(P), N)
  expect_equal(nrow(P), nrow(A_c))
  
  # Check that rows of P are not degenerate
  row_norms <- sqrt(rowSums(P^2))
  expect_true(all(row_norms > 1e-15))
  
  # Spectral check (loose)
  d <- rowSums(W)
  L <- Diagonal(x=d) - W
  orig_eigs <- sort(eigen(as.matrix(L), symmetric=TRUE, only.values=TRUE)$values)[1:2]
  
  d_c <- rowSums(A_c)
  L_c <- Diagonal(x=d_c) - A_c
  coars_eigs <- sort(eigen(as.matrix(L_c), symmetric=TRUE, only.values=TRUE)$values)[1:2]
  
  # Just check no pathological eigenvalues
  expect_true(all(is.finite(coars_eigs)))
  # Ensure no large negative values
  expect_true(all(coars_eigs >= -1e-12))
})

test_that("Coarsening does not produce invalid P or empty graphs", {
  # Use a slightly denser random graph
  set.seed(123)
  N <- 6
  A <- rsparsematrix(N, N, density=0.5, symmetric=TRUE)
  A <- (A + t(A))/2
  A[A<0] <- 0
  diag(A) <- 0

  epsilon_prime <- .5
  
  # Coarsen aiming at half size
  result <- coarsen_graph_multilevel(A, n=3, k=2, epsilon_prime=epsilon_prime, verbose=TRUE)
  
  A_c <- result$A_coarse
  P <- result$P
  info <- result$info
  
  # Basic checks
  expect_lt(nrow(A_c), N)
  expect_true(isSymmetric(A_c))
  expect_true(all(A_c@x >= 0))
  
  # P checks
  expect_equal(ncol(P), N)
  expect_equal(nrow(P), nrow(A_c))
  expect_true(all(rowSums(P) > 0)) # Each coarse vertex should have some contribution from original vertices
  
  # L_c computation
  d_c <- rowSums(A_c)
  L_c <- Diagonal(x=d_c) - A_c
  
  # Ensure not empty
  expect_gt(length(A_c@x), 0)
  
  # Eigenvalues finite
  coars_vals <- eigen(as.matrix(L_c), symmetric=TRUE, only.values=TRUE)$values
  expect_true(all(is.finite(coars_vals)))
  
  # Just a sanity check that epsilon and levels make sense
  # info$final_epsilon should be <= epsilon_prime or reflect that we stopped early
  expect_true(info$final_epsilon <= 1+epsilon_prime) # It's possible we ended early, but not expected to explode
  
  # Ensure that we actually did at least one contraction if possible
  expect_lt(nrow(A_c), N)
})


test_that("Coarsening a C_6 cycle graph down to 3 nodes yields a triangle with known spectrum", {
  # Create a cycle graph C_6
  # Vertices: 1 through 6
  # Edges: (1-2), (2-3), (3-4), (4-5), (5-6), (6-1)
  N <- 6
  A <- Matrix::sparseMatrix(
    i = c(1,2,3,4,5,6),
    j = c(2,3,4,5,6,1),
    x = rep(1, 6),
    dims = c(N,N),
    symmetric = TRUE
  )
  
  # Force to dgCMatrix class for consistency
  A <- as(A, "dgCMatrix")
  
  # We want to coarsen C_6 to n=3
  n_target <- 3
  k <- 2
  epsilon_prime <- 0.5
  
  result <- coarsen_graph_multilevel(A, n=n_target, k=k, epsilon_prime=epsilon_prime, verbose=FALSE)
  
  A_c <- result$A_coarse
  P <- result$P
  info <- result$info

  # Check final dimensions
  expect_equal(nrow(A_c), n_target)
  expect_equal(ncol(P), N)
  expect_equal(nrow(P), n_target)

  # Check that the coarsened adjacency is a triangle
  # We expect each of the 3 vertices to be connected to the other two with weight 1
  # and no self-loops.
  A_c_full <- as.matrix(A_c)
  expect_equal(dim(A_c_full), c(3,3))
  # Diagonal should be zero
  expect_true(all(abs(diag(A_c_full)) < 1e-14))
  # Off-diagonal should be 1
  off_diag <- A_c_full[row(A_c_full)!=col(A_c_full)]
  expect_true(all(abs(off_diag - 1) < 1e-14))

  # Compute the Laplacian of the coarsened graph
  d_c <- rowSums(A_c_full)
  L_c <- Diagonal(x=d_c) - A_c
  
  # Compute eigenvalues
  eigvals <- sort(eigs_sym(L_c, k=3, which="SM")$values)
  
  # Known eigenvalues for a triangle (K3/C3) are 0, 3, and 3
  expected_vals <- c(0,3,3)
  # Check within a tolerance
  expect_true(max(abs(eigvals - expected_vals)) < 1e-10)
})