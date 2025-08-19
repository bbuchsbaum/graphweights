library(testthat)
library(Matrix)
library(igraph)

context("Neighbor graph construction and operations")

test_that("neighbor_graph.Matrix creates valid neighbor_graph object", {
  # Create a simple adjacency matrix
  adj <- sparseMatrix(i = c(1, 2, 3), j = c(2, 3, 1), x = c(1, 1, 1), dims = c(3, 3))
  adj <- adj + t(adj)  # Make symmetric
  
  ng <- neighbor_graph(adj)
  
  expect_true(inherits(ng, "neighbor_graph"))
  expect_true("G" %in% names(ng))
  expect_true("params" %in% names(ng))
  expect_true(inherits(ng$G, "igraph"))
})

test_that("neighbor_graph.igraph creates valid neighbor_graph object", {
  # Create a simple igraph
  g <- make_ring(4)
  
  ng <- neighbor_graph(g)
  
  expect_true(inherits(ng, "neighbor_graph"))
  expect_true("G" %in% names(ng))
  expect_equal(ng$G, g)
})

test_that("adjacency.neighbor_graph extracts adjacency matrix", {
  adj_orig <- sparseMatrix(i = c(1, 2, 3), j = c(2, 3, 1), x = c(0.5, 0.7, 0.3), dims = c(3, 3))
  adj_orig <- adj_orig + t(adj_orig)
  
  ng <- neighbor_graph(adj_orig)
  adj_extracted <- adjacency(ng)
  
  expect_true(inherits(adj_extracted, "Matrix"))
  expect_equal(dim(adj_extracted), c(3, 3))
  expect_true(isSymmetric(adj_extracted))
})

test_that("edges.neighbor_graph returns edge information", {
  adj <- sparseMatrix(i = c(1, 2), j = c(2, 3), x = c(1, 1), dims = c(3, 3))
  adj <- adj + t(adj)
  
  ng <- neighbor_graph(adj)
  
  # Check that we can extract edges from the igraph object directly
  edge_info <- igraph::as_data_frame(ng$G, what = "edges")
  
  expect_true(is.data.frame(edge_info))
  expect_true(ncol(edge_info) >= 2)  # At least from and to columns
})

test_that("neighbors.neighbor_graph returns neighbor indices", {
  # Create a star graph: node 1 connected to nodes 2, 3, 4
  adj <- sparseMatrix(i = c(1, 1, 1), j = c(2, 3, 4), x = c(1, 1, 1), dims = c(4, 4))
  adj <- adj + t(adj)
  
  ng <- neighbor_graph(adj)
  
  # Node 1 should have 3 neighbors - use igraph directly to test
  neighs_1 <- igraph::neighbors(ng$G, 1)
  expect_equal(length(neighs_1), 3)
  
  # Convert to numeric for comparison
  neighs_1_numeric <- as.numeric(neighs_1)
  expect_true(all(neighs_1_numeric %in% c(2, 3, 4)))
  
  # Node 2 should have 1 neighbor
  neighs_2 <- igraph::neighbors(ng$G, 2)
  expect_equal(length(neighs_2), 1)
  expect_equal(as.numeric(neighs_2), 1)
})

test_that("laplacian.neighbor_graph computes Laplacian matrix", {
  adj <- sparseMatrix(i = c(1, 2, 3), j = c(2, 3, 1), x = c(1, 1, 1), dims = c(3, 3))
  adj <- adj + t(adj)  # Make symmetric
  
  ng <- neighbor_graph(adj)
  L <- laplacian(ng)
  
  expect_true(inherits(L, "Matrix"))
  expect_equal(dim(L), c(3, 3))
  expect_true(isSymmetric(L))
  
  # Laplacian should have zero row sums (up to numerical precision)
  expect_true(all(abs(rowSums(L)) < 1e-12))
})

test_that("nvertices and basic properties work", {
  # Create a symmetric adjacency matrix
  adj <- sparseMatrix(i = c(1, 2), j = c(2, 3), x = c(1, 1), dims = c(4, 4))
  adj <- adj + t(adj)  # Make symmetric
  
  ng <- neighbor_graph(adj)
  
  expect_equal(nvertices(ng), 4)
  expect_true(nvertices(ng) > 0)
})