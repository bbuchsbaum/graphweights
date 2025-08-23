library(testthat)
library(Matrix)

context("Diffusion and kernel methods")

test_that("heat_kernel produces valid weights", {
  distances <- c(0, 1, 2, 3, 4)
  sigma <- 1
  
  weights <- heat_kernel(distances, sigma = sigma)
  
  expect_equal(length(weights), 5)
  expect_true(all(weights >= 0))
  expect_true(all(weights <= 1))
  
  # Heat kernel should decrease with distance
  expect_true(weights[1] >= weights[2])
  expect_true(weights[2] >= weights[3])
  
  # At distance 0, weight should be 1
  expect_equal(weights[1], 1)
})

test_that("inverse_heat_kernel produces valid weights", {
  distances <- c(0, 1, 2, 3)
  sigma <- 1
  
  weights <- inverse_heat_kernel(distances, sigma = sigma)
  
  expect_equal(length(weights), 4)
  expect_true(all(weights >= 0))
  
  # Should handle zero distance appropriately
  expect_true(is.finite(weights[1]))
})

test_that("normalized_heat_kernel works correctly", {
  distances <- c(0, 1, 2)
  
  weights <- normalized_heat_kernel(distances, sigma = 1, len = 3)
  
  expect_equal(length(weights), 3)
  expect_true(all(weights >= 0))
  expect_true(all(is.finite(weights)))
})

test_that("dist_to_sim works with Matrix objects", {
  # Create a distance matrix
  dist_mat <- Matrix(c(0, 1, 2, 1, 0, 3, 2, 3, 0), nrow = 3, sparse = TRUE)
  
  # Test with heat method
  sim_mat <- dist_to_sim(dist_mat, method = "heat", sigma = 1)
  
  expect_true(inherits(sim_mat, "Matrix"))
  expect_equal(dim(sim_mat), c(3, 3))
  expect_true(all(sim_mat@x >= 0))
})

test_that("normalize_adjacency produces valid normalized matrix", {
  # Create a simple adjacency matrix
  adj <- sparseMatrix(i = c(1, 2, 3, 1), j = c(2, 3, 1, 3), x = c(2, 3, 1, 4), dims = c(3, 3))
  
  norm_adj <- normalize_adjacency(adj)
  
  expect_true(inherits(norm_adj, "Matrix"))
  expect_equal(dim(norm_adj), c(3, 3))
  expect_true(all(norm_adj@x >= 0))
  # Note: normalization might not bound values to [0,1] depending on method
  expect_true(all(is.finite(norm_adj@x)))
})

test_that("make_doubly_stochastic produces valid stochastic matrix", {
  # Create a symmetric positive matrix for better convergence
  adj <- sparseMatrix(i = c(1, 2, 3, 2, 3, 1), j = c(2, 3, 1, 1, 2, 3), 
                     x = c(1, 1, 1, 1, 1, 1), dims = c(3, 3))
  
  stoch_adj <- make_doubly_stochastic(adj)
  
  expect_true(inherits(stoch_adj, "Matrix"))
  expect_equal(dim(stoch_adj), c(3, 3))
  
  # Row sums should be approximately 1 (with reasonable tolerance)
  row_sums <- rowSums(stoch_adj)
  expect_true(all(abs(row_sums - 1) < 1e-6))
  
  # Column sums should be approximately 1 (with reasonable tolerance)  
  col_sums <- Matrix::colSums(stoch_adj)
  expect_true(all(abs(col_sums - 1) < 1e-6))
})

test_that("threshold_adjacency keeps top k values per row", {
  # Create adjacency matrix with various weights
  adj <- sparseMatrix(i = c(1, 1, 1, 2, 2, 3), j = c(1, 2, 3, 2, 3, 3), 
                     x = c(0.1, 0.5, 0.8, 0.3, 0.9, 0.2), dims = c(3, 3))
  
  # Keep top 2 values per row
  thresh_adj <- threshold_adjacency(adj, k = 2)
  
  expect_true(inherits(thresh_adj, "Matrix"))
  expect_equal(dim(thresh_adj), c(3, 3))
  
  # Each row should have at most k non-zero elements
  row_nnz <- Matrix::rowSums(thresh_adj != 0)
  expect_true(all(row_nnz <= 2))
})