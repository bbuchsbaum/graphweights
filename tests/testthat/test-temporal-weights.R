library(testthat)
library(Matrix)

context("Temporal weights and adjacency functions")

test_that("temporal_adjacency creates valid temporal connections", {
  # Simple temporal coordinates (time points)
  times <- c(1, 2, 3, 4, 5)
  
  result <- temporal_adjacency(times, window = 2)
  
  expect_true(inherits(result, "Matrix"))
  expect_equal(nrow(result), 5)
  expect_equal(ncol(result), 5)
  expect_true(all(result@x >= 0))
  
  # Convert to matrix for easier testing
  mat <- as.matrix(result)
  
  # Check that adjacent time points are connected
  expect_true(mat[1, 2] > 0)  # t=1 to t=2
  expect_true(mat[2, 3] > 0)  # t=2 to t=3
  
  # Check that non-adjacent time points are not connected (with lag=1)
  expect_equal(mat[1, 3], 0)  # t=1 to t=3 (lag > 1)
})

test_that("temporal_laplacian produces valid Laplacian", {
  times <- c(1, 2, 3, 4)
  
  result <- temporal_laplacian(times, window = 2)
  
  expect_true(inherits(result, "Matrix"))
  expect_equal(dim(result), c(4, 4))
  
  # Laplacian should have zero row sums
  row_sums <- rowSums(result)
  expect_true(all(abs(row_sums) < 1e-12))
})

test_that("temporal_autocor computes valid autocorrelation", {
  # Create simple time series data
  X <- matrix(1:20, nrow = 4, ncol = 5)  # 4 time points, 5 variables
  times <- c(1, 2, 3, 4)
  
  result <- temporal_autocor(X, window = 3)
  
  expect_true(inherits(result, "Matrix"))
  expect_equal(dim(result), c(5, 5))  # Correlates the 5 columns of X
  expect_true(all(result@x >= -1 & result@x <= 1))  # Correlation bounds
})

test_that("temporal functions handle edge cases", {
  # Single time point
  single_time <- c(1)
  result_single <- temporal_adjacency(single_time, window = 2)
  expect_equal(dim(result_single), c(1, 1))
  
  # Two time points
  two_times <- c(1, 2)
  result_two <- temporal_adjacency(two_times, window = 2)
  expect_equal(dim(result_two), c(2, 2))
  
  # Check that they are connected
  mat_two <- as.matrix(result_two)
  expect_true(mat_two[1, 2] > 0)
  expect_true(mat_two[2, 1] > 0)
})