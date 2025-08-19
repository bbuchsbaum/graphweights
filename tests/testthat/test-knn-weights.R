library(testthat)
library(Matrix)

context("KNN weights and neighbor functions")

test_that("weighted_knn produces valid adjacency matrix", {
  # Simple 2D feature matrix
  X <- matrix(c(1, 1,
                2, 2, 
                1, 2,
                2, 1,
                10, 10), ncol = 2, byrow = TRUE)
  
  # Test with sparse matrix output
  adj <- weighted_knn(X, k = 3, as = "sparse")
  
  # Basic structure tests
  expect_true(inherits(adj, "Matrix"))
  expect_equal(nrow(adj), 5)
  expect_equal(ncol(adj), 5)
  expect_true(all(adj@x >= 0))  # Non-negative weights
  
  # Should be sparse
  expect_true(inherits(adj, "sparseMatrix"))
})

test_that("weighted_knn respects k parameter", {
  X <- matrix(rnorm(20), ncol = 2)
  k <- 3
  
  adj <- weighted_knn(X, k = k, as = "sparse")
  
  # Each row should have approximately k non-zero elements (some variation allowed)
  row_nnz <- diff(adj@p)
  expect_true(all(row_nnz <= k + 2))  # Allow some flexibility
  expect_true(mean(row_nnz) <= k + 1)  # Average should be close to k
})

test_that("weighted_knn with different kernels", {
  X <- matrix(c(0, 0, 1, 0, 0, 1), ncol = 2, byrow = TRUE)
  
  # Test different kernel functions - don't pass sigma to rflann
  adj_heat <- weighted_knn(X, k = 2, FUN = heat_kernel, as = "sparse")
  adj_inv_heat <- weighted_knn(X, k = 2, FUN = inverse_heat_kernel, as = "sparse")
  
  expect_true(inherits(adj_heat, "Matrix"))
  expect_true(inherits(adj_inv_heat, "Matrix"))
  expect_equal(dim(adj_heat), c(3, 3))
  expect_equal(dim(adj_inv_heat), c(3, 3))
})

test_that("find_nn basic functionality", {
  X <- matrix(c(0, 0, 1, 0, 0, 1, 1, 1), ncol = 2, byrow = TRUE)
  
  # Create searcher object first - fix labels parameter
  searcher <- nnsearcher(X, labels = 1:nrow(X))
  
  # Find 2 nearest neighbors
  nn_result <- find_nn(searcher, k = 2)
  
  expect_true(is.list(nn_result))
  expect_true("indices" %in% names(nn_result))
  expect_true("distances" %in% names(nn_result))
  expect_equal(nrow(nn_result$indices), 4)  # 4 points
  expect_equal(nrow(nn_result$distances), 4)  # 4 points
})

test_that("find_nn_between works with different datasets", {
  X1 <- matrix(c(0, 0, 1, 0), ncol = 2, byrow = TRUE)
  X2 <- matrix(c(0, 1, 1, 1, 2, 1), ncol = 2, byrow = TRUE)
  
  # Create searcher object for X2 - fix labels parameter
  searcher <- nnsearcher(X2, labels = 1:nrow(X2))
  
  nn_result <- find_nn_between(searcher, k = 2, idx1 = 1:2, idx2 = 1:3)
  
  expect_true(is.list(nn_result))
  expect_true("indices" %in% names(nn_result))
  expect_true("distances" %in% names(nn_result))
})

test_that("heat_kernel and inverse_heat_kernel produce valid weights", {
  # Test distance-based weights directly
  distances <- c(0, 1, 2, 3)
  
  weights_heat <- heat_kernel(distances, sigma = 1)
  weights_inv <- inverse_heat_kernel(distances, sigma = 1)
  
  expect_equal(length(weights_heat), 4)
  expect_equal(length(weights_inv), 4)
  expect_true(all(weights_heat >= 0))
  expect_true(all(weights_inv >= 0))
  
  # Heat kernel should decrease with distance
  expect_true(weights_heat[1] >= weights_heat[2])
  expect_true(weights_heat[2] >= weights_heat[3])
})