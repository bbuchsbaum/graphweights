library(testthat)
library(Matrix)

context("Label similarity and adjacency functions")

test_that("binary_label_matrix works correctly for same labels", {
  labels <- c('A', 'A', 'B', 'B', 'C', 'A')
  
  result <- binary_label_matrix(labels, labels, type = 's')
  
  # Basic structure tests
  expect_true(inherits(result, "Matrix"))
  expect_equal(dim(result), c(6, 6))
  expect_true(all(result@x %in% c(0, 1)))  # Binary values
  
  # Convert to matrix for easier testing
  mat <- as.matrix(result)
  
  # Check that vertices with same labels are connected
  expect_equal(mat[1, 2], 1)  # A-A connection
  expect_equal(mat[1, 6], 1)  # A-A connection  
  expect_equal(mat[2, 6], 1)  # A-A connection
  expect_equal(mat[3, 4], 1)  # B-B connection
  
  # Check that vertices with different labels are NOT connected
  expect_equal(mat[1, 3], 0)  # A-B no connection
  expect_equal(mat[1, 5], 0)  # A-C no connection
  expect_equal(mat[3, 5], 0)  # B-C no connection
  
  # Check diagonal (self-connections)
  expect_equal(mat[1, 1], 1)  # Self-connection
  expect_equal(mat[5, 5], 1)  # Self-connection
})

test_that("binary_label_matrix works correctly for different labels", {
  labels <- c('A', 'A', 'B', 'B', 'C', 'A')
  
  result <- binary_label_matrix(labels, labels, type = 'd')
  
  # Basic structure tests
  expect_true(inherits(result, "Matrix"))
  expect_equal(dim(result), c(6, 6))
  
  # Convert to matrix for easier testing
  mat <- as.matrix(result)
  
  # Check that vertices with different labels are connected
  expect_equal(mat[1, 3], 1)  # A-B connection
  expect_equal(mat[1, 5], 1)  # A-C connection
  expect_equal(mat[3, 5], 1)  # B-C connection
  
  # Check that vertices with same labels are NOT connected
  expect_equal(mat[1, 2], 0)  # A-A no connection
  expect_equal(mat[1, 6], 0)  # A-A no connection
  expect_equal(mat[3, 4], 0)  # B-B no connection
  
  # Check diagonal (no self-connections for different type)
  expect_equal(mat[1, 1], 0)  # No self-connection
  expect_equal(mat[5, 5], 0)  # No self-connection
})

test_that("binary_label_matrix handles edge cases", {
  # Single label
  single_label <- c('A')
  result_single <- binary_label_matrix(single_label, single_label, type = 's')
  expect_equal(dim(result_single), c(1, 1))
  expect_equal(as.matrix(result_single)[1, 1], 1)
  
  # All different labels
  diff_labels <- c('A', 'B', 'C', 'D')
  result_diff <- binary_label_matrix(diff_labels, diff_labels, type = 's')
  mat_diff <- as.matrix(result_diff)
  # Should only have diagonal connections
  expect_equal(sum(mat_diff) - sum(diag(mat_diff)), 0)
  
  # All same labels
  same_labels <- c('A', 'A', 'A')
  result_same <- binary_label_matrix(same_labels, same_labels, type = 's')
  mat_same <- as.matrix(result_same)
  # Should be all 1s
  expect_true(all(mat_same == 1))
})

test_that("label_matrix works correctly", {
  labels1 <- c('A', 'B', 'A')
  labels2 <- c('A', 'A', 'B')
  
  result_same <- label_matrix(labels1, labels2, type = 's')
  result_diff <- label_matrix(labels1, labels2, type = 'd')
  
  expect_true(inherits(result_same, "Matrix"))
  expect_true(inherits(result_diff, "Matrix"))
  expect_equal(dim(result_same), c(3, 3))
  expect_equal(dim(result_diff), c(3, 3))
})

test_that("expand_label_similarity works with similarity matrix", {
  # Create a small similarity matrix
  sim_mat <- matrix(c(1.0, 0.8, 0.2,
                     0.8, 1.0, 0.3,
                     0.2, 0.3, 1.0), nrow = 3)
  rownames(sim_mat) <- colnames(sim_mat) <- c('A', 'B', 'C')
  
  labels <- c('A', 'B', 'A', 'C')
  
  # Test above threshold
  result_above <- expand_label_similarity(labels, sim_mat, threshold = 0.5, above = TRUE)
  expect_true(inherits(result_above, "Matrix"))
  expect_equal(dim(result_above), c(4, 4))
  
  # Test below threshold  
  result_below <- expand_label_similarity(labels, sim_mat, threshold = 0.5, above = FALSE)
  expect_true(inherits(result_below, "Matrix"))
  expect_equal(dim(result_below), c(4, 4))
})