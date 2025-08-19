library(testthat)
library(Matrix)

context("Repulsion and spatial constraint functions")

test_that("repulsion_graph creates valid repulsion matrix", {
  # Create simple coordinates and build spatial adjacency first
  coords <- matrix(c(0, 0, 1, 0, 0, 1, 1, 1), ncol = 2, byrow = TRUE)
  
  # Create a neighbor graph first
  W <- spatial_adjacency(coords, nnk = 3, sigma = 1, weight_mode = "heat")
  
  # Create class labels (alternating classes)
  labels <- factor(c(1, 2, 1, 2))
  cg <- class_graph(labels)
  
  result <- repulsion_graph(W, cg, method = "weighted")
  
  expect_true(inherits(result, "repulsion_graph"))
  expect_true(inherits(result, "neighbor_graph"))
  expect_equal(dim(adjacency(result)), c(4, 4))
  expect_true(all(adjacency(result)@x >= 0))  # Repulsion weights should be non-negative
})

test_that("spatial_constraints creates valid constraint matrix", {
  # Create coordinates
  coords <- matrix(c(0, 0, 1, 0, 0, 1, 1, 1), ncol = 2, byrow = TRUE)
  
  result <- spatial_constraints(coords, nblocks = 1, sigma_within = 1, nnk_within = 2)
  
  expect_true(inherits(result, "Matrix"))
  expect_equal(nrow(result), 4)
  expect_equal(ncol(result), 4)
})

test_that("feature_weighted_spatial_constraints works", {
  # Create coordinates and features for multiple blocks
  coords <- matrix(c(0, 0, 1, 0, 0, 1, 1, 1), ncol = 2, byrow = TRUE)
  features <- list(
    matrix(c(1, 1, 2, 2), nrow = 2, ncol = 4, byrow = TRUE),
    matrix(c(2, 2, 1, 1), nrow = 2, ncol = 4, byrow = TRUE)
  )
  
  result <- feature_weighted_spatial_constraints(coords, features, 
                                               sigma_within = 1, alpha_within = 0.5, nnk_within = 4)
  
  expect_true(inherits(result, "Matrix"))
  expect_equal(dim(result), c(8, 8))  # 4 coords * 2 blocks = 8
  expect_true(all(result@x >= 0))
})

test_that("spatial_smoother creates valid smoothing matrix", {
  # Create simple spatial coordinates
  coords <- matrix(c(0, 0, 1, 0, 0, 1), ncol = 2, byrow = TRUE)
  
  result <- spatial_smoother(coords, sigma = 1, nnk = 2)
  
  expect_true(inherits(result, "Matrix"))
  expect_equal(dim(result), c(3, 3))
  
  # Smoother should be row-stochastic (rows sum to 1)
  row_sums <- rowSums(result)
  expect_true(all(abs(row_sums - 1) < 1e-10))
})

test_that("bilateral_smoother works correctly", {
  # Create coordinates and features
  coords <- matrix(c(0, 0, 1, 0, 0, 1), ncol = 2, byrow = TRUE)
  features <- matrix(c(1, 1, 2, 2, 1, 2), ncol = 2, byrow = TRUE)
  
  result <- bilateral_smoother(coords, features, s_sigma = 1, f_sigma = 1, nnk = 4)
  
  expect_true(inherits(result, "Matrix"))
  expect_equal(dim(result), c(3, 3))
  
  # Should produce valid numeric results
  expect_true(all(is.finite(result@x)))
})

test_that("spatial_laplacian produces valid Laplacian", {
  coords <- matrix(c(0, 0, 1, 0, 0, 1, 1, 1), ncol = 2, byrow = TRUE)
  
  result <- spatial_laplacian(coords, dthresh = 2, nnk = 3, weight_mode = "heat", sigma = 1)
  
  expect_true(inherits(result, "Matrix"))
  expect_equal(dim(result), c(4, 4))
  expect_true(isSymmetric(result))
  
  # Laplacian should have zero row sums
  row_sums <- rowSums(result)
  expect_true(all(abs(row_sums) < 1e-12))
})

test_that("difference_of_gauss creates valid DoG filter", {
  coords <- matrix(c(0, 0, 1, 0, 0, 1, 1, 1), ncol = 2, byrow = TRUE)
  
  result <- difference_of_gauss(coords, sigma1 = 0.5, sigma2 = 1.5, nnk = 3)
  
  expect_true(inherits(result, "Matrix"))
  expect_equal(dim(result), c(4, 4))
})

test_that("spatial_lap_of_gauss creates valid LoG filter", {
  coords <- matrix(c(0, 0, 1, 0, 0, 1, 1, 1), ncol = 2, byrow = TRUE)
  
  result <- spatial_lap_of_gauss(coords, sigma = 1)
  
  expect_true(inherits(result, "Matrix"))
  expect_equal(dim(result), c(4, 4))
})