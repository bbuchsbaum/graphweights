test_that("class_graph constructor works correctly", {
  # Test with basic example
  labels <- factor(c("A", "A", "B", "B", "C"))
  cg <- class_graph(labels)
  
  # Check structure
  expect_type(cg, "list")
  expect_true("class_graph" %in% cg$classes)
  expect_equal(length(cg$levels), 3)
  expect_equal(cg$levels, c("A", "B", "C"))
  
  # Check class frequencies
  expect_equal(as.vector(cg$class_freq), c(2, 2, 1))
  
  # Check class indices
  expect_equal(cg$class_indices$A, c(1, 2))
  expect_equal(cg$class_indices$B, c(3, 4))
  expect_equal(cg$class_indices$C, 5)
  
  # Check adjacency matrix structure
  expect_true(methods::is(cg$adjacency, "Matrix"))
  expect_equal(dim(cg$adjacency), c(5, 5))
  
  # Check that same-class entries are 1, different-class are 0
  expect_equal(as.vector(cg$adjacency[1,2]), 1)  # A-A connection
  expect_equal(as.vector(cg$adjacency[3,4]), 1)  # B-B connection
  expect_equal(as.vector(cg$adjacency[1,3]), 0)  # A-B connection
})

test_that("nclasses.class_graph works correctly", {
  labels <- factor(c("A", "A", "B", "B", "C"))
  cg <- class_graph(labels)
  expect_equal(nclasses.class_graph(cg), 3)
  
  # Test with single class
  cg_single <- class_graph(factor(rep("A", 3)))
  expect_equal(nclasses.class_graph(cg_single), 1)
})

test_that("class_means.class_graph works correctly", {
  # Create test data
  X <- matrix(1:12, nrow=4, ncol=3)
  labels <- factor(c("A", "A", "B", "B"))
  cg <- class_graph(labels)
  
  means <- class_means.class_graph(cg, X)
  
  # Check structure
  expect_equal(dim(means), c(2, 3))
  expect_equal(rownames(means), c("A", "B"))
  
  # Check actual means
  expect_equal(means["A",], colMeans(X[1:2,]))
  expect_equal(means["B",], colMeans(X[3:4,]))
})

test_that("within_class_neighbors.class_graph works correctly", {
  # Create test data
  X <- matrix(1:12, nrow=4, ncol=3)
  labels <- factor(c("A", "A", "B", "B"))
  cg <- class_graph(labels)
  
  # Create neighbor graph using weighted_knn
  adj <- weighted_knn(X, k=2, as="sparse")
  ng <- neighbor_graph(adj)
  
  
  within_ng <- within_class_neighbors.class_graph(cg, ng)
  
  # Check structure
  expect_true(methods::is(adjacency(within_ng), "Matrix"))
  expect_equal(dim(adjacency(within_ng)), c(4, 4))
  
  # Check that only within-class connections exist
  adj <- as.matrix(adjacency(within_ng))
  expect_true(all(adj[1:2, 3:4] == 0))  # No A-B connections
})

test_that("between_class_neighbors.class_graph works correctly", {
  # Create test data
  X <- matrix(1:12, nrow=4, ncol=3)
  labels <- factor(c("A", "A", "B", "B"))
  cg <- class_graph(labels)
  
  # Create neighbor graph using weighted_knn
  adj <- weighted_knn(X, k=2, as="sparse")
  ng <- neighbor_graph(adj)
  
  between_ng <- between_class_neighbors.class_graph(cg, ng)
  
  # Check structure
  expect_true(methods::is(adjacency(between_ng), "Matrix"))
  expect_equal(dim(adjacency(between_ng)), c(4, 4))
  
  # Check that only between-class connections exist
  adj <- as.matrix(adjacency(between_ng))
  expect_true(all(adj[1:2, 1:2] == 0))  # No within-A connections
  expect_true(all(adj[3:4, 3:4] == 0))  # No within-B connections
})

test_that("discriminating_distance works correctly", {
  # Create test data
  set.seed(123)
  X <- matrix(rnorm(20), nrow=4, ncol=5)
  labels <- factor(c("A", "A", "B", "B"))
  
  D <- discriminating_distance(X, k=2, sigma=0.5, labels=labels)
  
  # Check structure
  expect_true(is.matrix(D))
  expect_equal(dim(D), c(4, 4))
  
  # Check properties
  expect_true(isSymmetric(D))  # Should be symmetric
  expect_true(all(diag(D) == 0))  # Diagonal should be zero
  expect_true(all(D >= 0))  # Distances should be non-negative
})

test_that("class_graph handles edge cases", {
  # Single observation
  expect_error(class_graph(factor("A")), NA)
  
  # Empty labels
  expect_error(class_graph(factor(character(0))))
  
  # NA labels
  expect_error(class_graph(factor(c("A", NA, "B"))))
  
  # Non-factor input
  labels <- c("A", "B", "C")
  cg <- class_graph(labels)
  expect_true(is.factor(cg$labels))
})

test_that("class_graph handles self-connections correctly", {
  labels <- factor(c("A", "A", "B", "B", "C"))
  
  # Test with self_connections = FALSE (default)
  cg_no_self <- class_graph(labels)
  expect_equal(diag(cg_no_self$adjacency), rep(0, 5))
  
  # Test with self_connections = TRUE
  cg_with_self <- class_graph(labels, self_connections = TRUE)
  expect_equal(diag(cg_with_self$adjacency), rep(1, 5))
  
  # Verify other connections remain the same
  adj_no_self <- as.matrix(cg_no_self$adjacency)
  adj_with_self <- as.matrix(cg_with_self$adjacency)
  
  # Zero out diagonals for comparison
  diag(adj_no_self) <- 0
  diag(adj_with_self) <- 0
  
  expect_equal(adj_no_self, adj_with_self)
})

test_that("print.class_graph works correctly", {
  # Create test data
  labels <- factor(c("A", "A", "B", "B", "C"))
  cg <- class_graph(labels)
  
  # Capture the output
  output <- capture.output(print(cg))
  
  # Basic checks
  expect_true(any(grepl("Class Graph Summary", output)))
  expect_true(any(grepl("Basic Information", output)))
  expect_true(any(grepl("Class Distribution", output)))
  expect_true(any(grepl("Edge Information", output)))
  
  # Check content
  expect_true(any(grepl("5", output)))  # Number of vertices
  expect_true(any(grepl("3", output)))  # Number of classes
  expect_true(any(grepl("A", output)))  # Class names
  expect_true(any(grepl("B", output)))
  expect_true(any(grepl("C", output)))
})