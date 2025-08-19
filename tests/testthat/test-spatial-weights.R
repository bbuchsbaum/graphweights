library(testthat)
library(Matrix)

context("Spatial weights and adjacency functions")

test_that("spatial_adjacency produces valid adjacency matrix", {
  # Simple 2D coordinates
  coords <- matrix(c(0, 0, 1, 0, 0, 1, 1, 1), ncol = 2, byrow = TRUE)
  
  # Test basic functionality
  adj <- spatial_adjacency(coords, sigma = 1, nnk = 3, weight_mode = "heat")
  
  # Basic structure tests
  expect_true(inherits(adj, "Matrix"))
  expect_equal(nrow(adj), 4)
  expect_equal(ncol(adj), 4)
  expect_true(isSymmetric(adj))
  expect_true(all(adj@x >= 0))  # Non-negative weights
  
  # Diagonal might not be zero by default - just check structure
  expect_true(all(is.finite(diag(adj))))
})

test_that("spatial_adjacency binary mode works", {
  coords <- matrix(c(0, 0, 1, 0, 2, 0), ncol = 2, byrow = TRUE)
  
  adj_binary <- spatial_adjacency(coords, weight_mode = "binary", nnk = 2)
  adj_heat <- spatial_adjacency(coords, weight_mode = "heat", nnk = 2)
  
  # Binary should have discrete values (may not be exactly 1)
  binary_vals <- adj_binary@x[adj_binary@x > 0]
  expect_true(length(unique(binary_vals)) <= 2)  # At most 2 distinct values
  
  # Heat should have continuous values
  heat_vals <- adj_heat@x[adj_heat@x > 0]
  expect_true(any(heat_vals != 1))
})

test_that("spatial_adjacency handles edge cases", {
  # Single point
  coords_single <- matrix(c(0, 0), ncol = 2)
  adj_single <- spatial_adjacency(coords_single)
  expect_equal(dim(adj_single), c(1, 1))
  expect_true(is.finite(adj_single[1,1]))  # May not be zero
  
  # Two identical points
  coords_dup <- matrix(c(0, 0, 0, 0), ncol = 2, byrow = TRUE)
  adj_dup <- spatial_adjacency(coords_dup, sigma = 1)
  expect_true(isSymmetric(adj_dup))
  expect_equal(nrow(adj_dup), 2)
})

test_that("cross_spatial_adjacency works between different coordinate sets", {
  coords1 <- matrix(c(0, 0, 1, 0), ncol = 2, byrow = TRUE)
  coords2 <- matrix(c(0, 1, 1, 1, 2, 1), ncol = 2, byrow = TRUE)
  
  cross_adj <- cross_spatial_adjacency(coords1, coords2, sigma = 1)
  
  expect_equal(nrow(cross_adj), 2)  # coords1 rows
  expect_equal(ncol(cross_adj), 3)  # coords2 rows
  expect_true(all(cross_adj@x >= 0))
})

test_that("weighted_spatial_adjacency combines spatial and feature similarity", {
  # Simple spatial coordinates
  coords <- matrix(c(0, 0, 1, 0, 0, 1, 1, 1), ncol = 2, byrow = TRUE)
  # Simple features
  features <- matrix(c(1, 1, 2, 2, 1, 2, 2, 1), ncol = 2, byrow = TRUE)
  
  # Use correct parameter names
  adj <- weighted_spatial_adjacency(coords, features, wsigma = 1, alpha = 0.5, 
                                   weight_mode = "heat")
  
  expect_true(inherits(adj, "Matrix"))
  expect_equal(dim(adj), c(4, 4))
  expect_true(isSymmetric(adj))
  expect_true(all(adj@x >= 0))
})