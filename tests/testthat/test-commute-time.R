library(testthat)
library(Matrix)

context("Commute time and diffusion functions")

test_that("commute_time_distance computes valid distances", {
  # Create a simple path graph: 1-2-3-4
  adj <- sparseMatrix(i = c(1, 2, 3), j = c(2, 3, 4), x = c(1, 1, 1), dims = c(4, 4))
  adj <- adj + t(adj)  # Make symmetric
  
  result <- commute_time_distance(adj)
  
  expect_true(inherits(result, "commute_time"))
  expect_true(is.list(result))
  expect_true("cds" %in% names(result))
  expect_true("eigenvalues" %in% names(result))
  expect_true("eigenvectors" %in% names(result))
  expect_true("gap" %in% names(result))
  
  # Check the cds matrix (embedding coordinates)
  expect_true(is.matrix(result$cds))
  expect_equal(nrow(result$cds), 4)  # Should match number of nodes
  
  # Check eigenvalues are valid
  expect_true(all(result$eigenvalues <= 1 + 1e-12))
  expect_true(is.numeric(result$gap))
})

test_that("heat_kernel produces valid similarity values", {
  distances <- c(0, 0.5, 1, 2, 5)
  sigma <- 1
  
  similarities <- heat_kernel(distances, sigma = sigma)
  
  expect_equal(length(similarities), 5)
  expect_true(all(similarities >= 0))
  expect_true(all(similarities <= 1))
  
  # Distance 0 should give similarity 1
  expect_equal(similarities[1], 1)
  
  # Similarities should decrease with distance
  expect_true(similarities[1] >= similarities[2])
  expect_true(similarities[2] >= similarities[3])
  expect_true(similarities[3] >= similarities[4])
})

test_that("inverse_heat_kernel handles zero distances", {
  distances <- c(0, 1, 2, 3)
  sigma <- 1
  
  similarities <- inverse_heat_kernel(distances, sigma = sigma)
  
  expect_equal(length(similarities), 4)
  expect_true(all(similarities >= 0))
  expect_true(all(is.finite(similarities)))
})

test_that("estimate_sigma provides reasonable estimates", {
  # Create a sample data matrix
  X <- matrix(rnorm(50), nrow = 10, ncol = 5)
  
  sigma_est <- estimate_sigma(X)
  
  expect_true(is.numeric(sigma_est))
  expect_true(length(sigma_est) == 1)
  expect_true(sigma_est > 0)
})

test_that("normalized_heat_kernel works with length parameter", {
  distances <- c(0, 1, 2, 3)
  sigma <- 1
  len <- 4
  
  similarities <- normalized_heat_kernel(distances, sigma = sigma, len = len)
  
  expect_equal(length(similarities), 4)
  expect_true(all(similarities >= 0))
  expect_true(all(is.finite(similarities)))
})