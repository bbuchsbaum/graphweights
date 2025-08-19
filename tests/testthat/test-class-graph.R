library(testthat)
library(Matrix)

context("Class graph construction and operations")

test_that("class_graph creates valid class-based graph", {
  # Create sample data with classes
  X <- matrix(rnorm(20), nrow = 10, ncol = 2)
  classes <- factor(c(1, 1, 2, 2, 3, 3, 1, 2, 3, 1))
  
  cg <- class_graph(classes)
  
  expect_true(inherits(cg, "class_graph"))
  expect_true("G" %in% names(cg))
  expect_true("labels" %in% names(cg))
  expect_equal(cg$labels, classes)
})

test_that("within_class_neighbors works correctly", {
  # Create sample data and simple adjacency matrix
  classes <- factor(c('A', 'A', 'B', 'B'))
  
  cg <- class_graph(classes)
  
  # Create a simple adjacency matrix (fully connected)
  library(Matrix)
  adj <- sparseMatrix(i = c(1, 1, 2, 2, 3, 3, 4, 4), 
                      j = c(2, 3, 1, 4, 1, 4, 2, 3), 
                      x = rep(1, 8), dims = c(4, 4))
  ng <- neighbor_graph(adj)
  
  # Find within-class neighbors
  within_neighs <- within_class_neighbors(cg, ng)
  
  expect_true(inherits(within_neighs, "neighbor_graph"))
  expect_true("G" %in% names(within_neighs))
})

test_that("between_class_neighbors works correctly", {
  # Create sample data and simple adjacency matrix
  classes <- factor(c('A', 'A', 'B', 'B'))
  
  cg <- class_graph(classes)
  
  # Create a simple adjacency matrix (fully connected)
  library(Matrix)
  adj <- sparseMatrix(i = c(1, 1, 2, 2, 3, 3, 4, 4), 
                      j = c(2, 3, 1, 4, 1, 4, 2, 3), 
                      x = rep(1, 8), dims = c(4, 4))
  ng <- neighbor_graph(adj)
  
  # Find between-class neighbors
  between_neighs <- between_class_neighbors(cg, ng)
  
  expect_true(inherits(between_neighs, "neighbor_graph"))
  expect_true("G" %in% names(between_neighs))
})

test_that("nclasses returns correct number of classes", {
  classes <- factor(c('A', 'A', 'B', 'B', 'C'))
  X <- matrix(rnorm(10), nrow = 5, ncol = 2)
  
  cg <- class_graph(classes)
  
  expect_equal(nclasses(cg), 3)
})

test_that("class_means computes correct centroids", {
  # Create simple data where means are obvious
  X <- matrix(c(1, 1,   # Class A point 1
                2, 2,   # Class A point 2  
                10, 10, # Class B point 1
                11, 11), # Class B point 2
              ncol = 2, byrow = TRUE)
  classes <- factor(c('A', 'A', 'B', 'B'))
  
  cg <- class_graph(classes)
  means <- class_means(cg, X)
  
  expect_true(is.matrix(means))
  expect_equal(nrow(means), 2)  # 2 classes
  expect_equal(ncol(means), 2)  # 2 dimensions
  
  # Check that means are approximately correct
  expect_true(abs(means[1, 1] - 1.5) < 0.1)  # Class A mean X
  expect_true(abs(means[1, 2] - 1.5) < 0.1)  # Class A mean Y
  expect_true(abs(means[2, 1] - 10.5) < 0.1) # Class B mean X
  expect_true(abs(means[2, 2] - 10.5) < 0.1) # Class B mean Y
})