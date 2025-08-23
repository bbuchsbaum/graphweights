library(testthat)
library(Matrix)

context("Design similarity kernels")

test_that("design_kernel creates valid kernel for nominal factors", {
  # Simple 2x3 design
  factors <- list(
    A = list(L=2, type="nominal"),
    B = list(L=3, type="nominal")
  )
  
  result <- design_kernel(factors)
  
  # Check structure
  expect_true(is.list(result))
  expect_true("K" %in% names(result))
  expect_true("K_cell" %in% names(result))
  expect_true("info" %in% names(result))
  
  # Check dimensions (2x3 = 6 cells)
  expect_equal(nrow(result$K), 6)
  expect_equal(ncol(result$K), 6)
  expect_equal(dim(result$K_cell), c(6, 6))
  
  # Check kernel properties
  K <- result$K_cell
  expect_true(Matrix::isSymmetric(K))
  expect_true(all(K >= -1e-10))  # Non-negative (allowing small numerical errors)
  
  # Check positive semi-definite
  eigs <- eigen(K, symmetric=TRUE, only.values=TRUE)$values
  expect_true(all(eigs >= -1e-10))
  
  # Check info structure
  expect_equal(result$info$levels, c(A=2, B=3))
  expect_equal(result$info$factor_names, c("A", "B"))
  expect_equal(result$info$basis, "cell")
})

test_that("design_kernel handles ordinal factors correctly", {
  factors <- list(
    dose = list(L=5, type="ordinal", l=1.5),
    group = list(L=2, type="nominal")
  )
  
  result <- design_kernel(factors)
  K <- result$K_cell
  
  # Check dimensions (5x2 = 10 cells)
  expect_equal(dim(K), c(10, 10))
  
  # Check that ordinal factor creates smooth transitions
  # Extract the dose kernel pattern (first 5x5 block structure should show smoothness)
  # Due to Kronecker structure, cells with adjacent dose levels should be more similar
  # This is a basic sanity check
  expect_true(Matrix::isSymmetric(K))
  expect_true(all(eigen(K, symmetric=TRUE, only.values=TRUE)$values >= -1e-10))
})

test_that("design_kernel handles circular factors correctly", {
  factors <- list(
    hour = list(L=4, type="circular", l=1.0)
  )
  
  result <- design_kernel(factors)
  K <- result$K
  
  # Check wrap-around: distance from level 1 to level 4 should equal distance from 1 to 2
  # Due to circular topology
  expect_equal(dim(K), c(4, 4))
  expect_true(Matrix::isSymmetric(K))
  
  # In a circular kernel, K[1,4] should be similar to K[1,2] due to wrap-around
  # The exact values depend on length scale, but symmetry should hold
  expect_true(abs(K[1,4] - K[2,3]) < 1e-10)  # Circular symmetry
})

test_that("design_kernel with custom terms works", {
  factors <- list(
    A = list(L=2, type="nominal"),
    B = list(L=2, type="nominal")
  )
  
  # Only main effects
  terms <- list("A", "B")
  rho <- c(A=1, B=1)
  
  result <- design_kernel(factors, terms=terms, rho=rho, include_intercept=FALSE)
  
  expect_equal(length(result$info$term_names), 2)
  expect_equal(result$info$term_names, c("A", "B"))
})

test_that("design_kernel with interaction terms works", {
  factors <- list(
    A = list(L=2, type="nominal"),
    B = list(L=2, type="nominal")
  )
  
  # Main effects and interaction
  terms <- list("A", "B", c("A", "B"))
  rho <- c(A=1, B=1, "A:B"=0.5)
  
  result <- design_kernel(factors, terms=terms, rho=rho)
  
  expect_equal(length(result$info$term_names), 3)
  expect_true("A:B" %in% result$info$term_names)
})

test_that("design_kernel effect-coding transformation works", {
  factors <- list(
    A = list(L=3, type="nominal"),
    B = list(L=2, type="nominal")
  )
  
  # Create contrasts
  contrasts <- list(
    A = contr.sum(3),
    B = contr.sum(2)
  )
  
  result <- design_kernel(factors, basis="effect", contrasts=contrasts)
  
  # Check that we're in effect space
  expect_equal(result$info$basis, "effect")
  
  # Effect space dimension should be smaller than cell space
  # A has 3-1=2 contrasts, B has 2-1=1 contrast
  # So we expect dimensions based on the terms included
  expect_true(!is.null(result$info$map))
  expect_true(!is.null(result$info$blocks))
  
  # K should still be PSD
  eigs <- eigen(result$K, symmetric=TRUE, only.values=TRUE)$values
  expect_true(all(eigs >= -1e-10))
})

test_that("kernel normalization options work", {
  factors <- list(A = list(L=3, type="nominal"))
  
  # Test unit trace normalization
  result_trace <- design_kernel(factors, normalize="unit_trace")
  expect_equal(sum(diag(result_trace$K)), 1, tolerance=1e-10)
  
  # Test max diagonal normalization
  result_maxdiag <- design_kernel(factors, normalize="max_diag")
  expect_equal(max(diag(result_maxdiag$K)), 1, tolerance=1e-10)
  
  # Test Frobenius normalization
  result_fro <- design_kernel(factors, normalize="unit_fro")
  fro_norm <- sqrt(sum(result_fro$K^2))
  expect_equal(fro_norm, 1, tolerance=1e-10)
})

test_that("kernel_roots computes valid square roots", {
  # Create a simple kernel
  factors <- list(A = list(L=3, type="nominal"))
  K_result <- design_kernel(factors)
  K <- K_result$K
  
  # Compute roots
  roots <- kernel_roots(K)
  
  expect_true(is.list(roots))
  expect_true(all(c("Khalf", "Kihalf", "evals", "evecs") %in% names(roots)))
  
  # Check that Khalf * Khalf = K (approximately)
  K_reconstructed <- roots$Khalf %*% roots$Khalf
  expect_equal(K, K_reconstructed, tolerance=1e-10)
  
  # Check that Kihalf * K * Kihalf = I (approximately)
  I_reconstructed <- roots$Kihalf %*% K %*% roots$Kihalf
  n <- nrow(K)
  expect_equal(I_reconstructed, diag(n), tolerance=1e-10)
})

test_that("kernel_alignment computes correct scores", {
  # Two identical matrices should have alignment 1
  A <- matrix(c(1,0.5,0.5,1), 2, 2)
  expect_equal(kernel_alignment(A, A), 1, tolerance=1e-10)
  
  # Orthogonal matrices should have alignment 0
  B <- matrix(c(1,0,0,1), 2, 2)
  C <- matrix(c(0,1,1,0), 2, 2)
  expect_equal(kernel_alignment(B, C), 0, tolerance=1e-10)
  
  # Check range [-1,1] (can be negative for negatively correlated matrices)
  D <- matrix(rnorm(16), 4, 4)
  D <- (D + t(D))/2  # Make symmetric
  E <- matrix(rnorm(16), 4, 4)
  E <- (E + t(E))/2
  
  alignment <- kernel_alignment(D, E)
  expect_true(alignment >= -1 - 1e-10 && alignment <= 1 + 1e-10)
})

test_that("example_kernel_5x5 creates correct structure", {
  # Test with default parameters
  result <- example_kernel_5x5()
  
  expect_true(is.list(result))
  expect_equal(dim(result$K), c(25, 25))
  expect_equal(dim(result$K_cell), c(25, 25))
  
  # Test with custom weights
  result2 <- example_kernel_5x5(rhoA=2, rhoB=3, rhoAB=1)
  expect_true(Matrix::isSymmetric(result2$K))
  
  # Check PSD
  eigs <- eigen(result2$K, symmetric=TRUE, only.values=TRUE)$values
  expect_true(all(eigs >= -1e-10))
})

test_that("sum_contrasts creates valid contrast matrices", {
  Ls <- c(A=3, B=4, C=2)
  contrasts <- sum_contrasts(Ls)
  
  expect_true(is.list(contrasts))
  expect_equal(names(contrasts), names(Ls))
  
  # Check dimensions
  expect_equal(dim(contrasts$A), c(3, 2))  # L x (L-1)
  expect_equal(dim(contrasts$B), c(4, 3))
  expect_equal(dim(contrasts$C), c(2, 1))
  
  # Check sum-to-zero property
  expect_equal(sum(contrasts$A[,1]), 0, tolerance=1e-10)
  expect_equal(sum(contrasts$B[,1]), 0, tolerance=1e-10)
})

test_that("helmert_contrasts creates orthonormal matrices", {
  Ls <- c(A=3, B=4)
  contrasts <- helmert_contrasts(Ls)
  
  expect_true(is.list(contrasts))
  expect_equal(names(contrasts), names(Ls))
  
  # Check dimensions
  expect_equal(dim(contrasts$A), c(3, 2))
  expect_equal(dim(contrasts$B), c(4, 3))
  
  # Check orthonormality
  A_gram <- t(contrasts$A) %*% contrasts$A
  expect_equal(A_gram, diag(2), tolerance=1e-10)
  
  B_gram <- t(contrasts$B) %*% contrasts$B
  expect_equal(B_gram, diag(3), tolerance=1e-10)
})

test_that("design_kernel handles edge cases", {
  # Single factor
  factors <- list(A = list(L=1, type="nominal"))
  result <- design_kernel(factors)
  expect_equal(dim(result$K), c(1, 1))
  
  # Large number of levels
  factors <- list(A = list(L=10, type="nominal"))
  result <- design_kernel(factors)
  expect_equal(dim(result$K), c(10, 10))
  expect_true(Matrix::isSymmetric(result$K))
})

test_that("design_kernel input validation works", {
  # Empty factors
  expect_error(design_kernel(list()), "non-empty")
  
  # Unnamed factors
  expect_error(design_kernel(list(list(L=2))), "NAMED list")
  
  # Missing levels
  expect_error(design_kernel(list(A=list(type="nominal"))), "needs L")
  
  # Invalid level count
  expect_error(design_kernel(list(A=list(L=0, type="nominal"))), "positive integers")
  
  # Invalid type
  expect_error(design_kernel(list(A=list(L=2, type="invalid"))), "nominal.*ordinal.*circular")
  
  # Duplicate factor names
  expect_error(design_kernel(list(A=list(L=2), A=list(L=3))), "unique")
  
  # Invalid length scale
  expect_error(design_kernel(list(A=list(L=3, type="ordinal", l=-1))), "positive")
  
  # Negative rho
  expect_error(design_kernel(list(A=list(L=2)), rho=c(A=-1)), "nonnegative")
  
  # Missing contrasts for effect basis
  expect_error(design_kernel(list(A=list(L=2)), basis="effect"), "requires.*contrasts")
})

test_that("design_kernel preserves deterministic output", {
  factors <- list(
    A = list(L=3, type="nominal"),
    B = list(L=2, type="ordinal", l=1.0)
  )
  
  # Run twice with same parameters
  result1 <- design_kernel(factors, terms=list("A", "B", c("A","B")), rho=c(A=1, B=2, "A:B"=0.5))
  result2 <- design_kernel(factors, terms=list("A", "B", c("A","B")), rho=c(A=1, B=2, "A:B"=0.5))
  
  # Results should be identical
  expect_equal(result1$K, result2$K)
  expect_equal(result1$K_cell, result2$K_cell)
  expect_equal(result1$info$term_names, result2$info$term_names)
})