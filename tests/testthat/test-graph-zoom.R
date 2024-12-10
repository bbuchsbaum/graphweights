# tests/testthat/test_graphzoom.R
library(testthat)
library(Matrix)

test_that("GraphZoom pipeline runs and produces correct output", {
  # Set seed for reproducibility
  set.seed(123)
  
  # Create a small synthetic graph
  N <- 30
  # Sparse random graph adjacency
  A <- rsparsematrix(N, N, density = 0.05)
  A <- (A + t(A)) / 2
  diag(A) <- 0
  
  # Create a small feature matrix
  X <- matrix(rnorm(N*5), nrow = N, ncol = 5)
  
  # Parameters
  levels <- 2
  embed_dim <- 4
  kNN <- 5
  filter_power <- 2
  sigma <- 1e-3
  
  # Run the GraphZoom pipeline
  # This should use the spectral embedding method at the coarsest graph if implemented.
  result <- graph_zoom(
    A = A,
    X = X,
    levels = levels,
    kNN = kNN,
    embed_dim = embed_dim,
    filter_power = filter_power,
    sigma = sigma,
    verbose = TRUE
  )
  
  # Checks
  # 1. Result contains expected fields
  expect_named(result, c("E_original", "graphs", "H_list", "E_coarse"))
  
  # 2. Check that graphs is a list of length levels+1
  # graphs[[1]] should be the fused/original-level graph
  # graphs[[levels+1]] is the coarsest graph
  expect_length(result$graphs, levels + 1)
  
  # 3. Check graph sizes decreasing
  sizes <- sapply(result$graphs, nrow)
  expect_true(all(diff(sizes) < 0))  # sizes should strictly decrease
  
  # 4. Check embedding dimensions at the end
  E_original <- result$E_original
  expect_equal(dim(E_original), c(N, embed_dim))
  
  # 5. Check embeddings are numeric and finite
  expect_true(is.numeric(as.matrix(E_original)))
  expect_true(all(is.finite(as.matrix(E_original))))
  
  # 6. Check coarsest embedding dimension
  E_coarse <- result$E_coarse
  # Coarsest graph size:
  coarsest_size <- nrow(result$graphs[[levels+1]])
  expect_equal(dim(E_coarse), c(coarsest_size, embed_dim))
  
  # 7. Check H_list length
  expect_length(result$H_list, levels)
  
  # This test primarily ensures that the pipeline runs without error,
  # that coarsening works as expected, and embeddings are produced 
  # with the correct size and no obvious anomalies.
  
  # If you want additional diagnostics, you could print summary stats,
  # but typically tests should remain silent unless there is a failure.
  
  # If spectral embedding was implemented, we trust the eigen-decomposition 
  # gives a meaningful embedding. We do not test for correctness of spectral 
  # embedding itself here, just that it runs and returns correct dimensions.
})

# tests/testthat/test_graphzoom_diagnostic.R

library(testthat)
library(Matrix)

test_that("GraphZoom diagnostic test on a simple, known graph", {
  # Construct a small, simple graph:
  # Let's use a simple chain of 5 nodes:
  # 1 -- 2 -- 3 -- 4 -- 5
  # Adjacency for a chain:
  N <- 5
  i <- c(1,2,3,4)
  j <- c(2,3,4,5)
  x <- rep(1, 4)
  A <- sparseMatrix(i=c(i,j), j=c(j,i), x=1, dims=c(N,N))
  
  # Node features: simple random vectors. For a small test, random is fine.
  # We'll add a small dimension of features, say 3:
  set.seed(42)
  X <- matrix(rnorm(N*3), nrow=N, ncol=3)
  
  # Parameters:
  levels <- 1  # Just one level of coarsening for simplicity
  kNN <- 2     # With a chain and small N, 2-NN from features is trivial
  embed_dim <- 2
  filter_power <- 2
  sigma <- 1e-3
  
  # Run the pipeline
  result <- graph_zoom(
    A = A,
    X = X,
    levels = levels,
    kNN = kNN,
    embed_dim = embed_dim,
    filter_power = filter_power,
    sigma = sigma,
    verbose = FALSE
  )
  
  # Check the structure of the result
  expect_named(result, c("E_original", "graphs", "H_list", "E_coarse"))
  
  # Check number of graphs: original + one coarsened
  expect_length(result$graphs, levels + 1)
  
  # The coarsened graph should have fewer nodes than original
  original_size <- nrow(result$graphs[[1]])
  coarsest_size <- nrow(result$graphs[[levels+1]])
  expect_true(coarsest_size < original_size)
  
  # Check embeddings dimension
  E_original <- result$E_original
  expect_equal(dim(E_original), c(N, embed_dim))
  
  # Check no NaNs or Infs in original embeddings
  expect_true(all(is.finite(E_original)))
  
  # Check coarsest embedding size and no NaNs
  E_coarse <- result$E_coarse
  expect_equal(dim(E_coarse), c(coarsest_size, embed_dim))
  expect_true(all(is.finite(E_coarse)))
  
  # Check that H_list is correct
  # There's 1 level, so 1 H:
  expect_length(result$H_list, levels)
  H <- result$H_list[[1]]
  # Dimensions: #coarse_nodes x #original_nodes
  expect_equal(dim(H), c(coarsest_size, original_size))
  
  # Optional: Since this is a chain graph, we expect that coarsening merges nodes.
  # Just a quick check that coarsest graph is smaller:
  expect_lt(coarsest_size, original_size)
  
  # The test ensures the pipeline runs and output is well-formed on a small, known graph.
  # This serves as a regression test if future changes break assumptions.
})


test_that("Providing A_feat or X yields the same result", {
  set.seed(42)
  
  # Create a small graph
  N <- 20
  A <- rsparsematrix(N, N, 0.05)
  A <- (A + t(A))/2
  diag(A) <- 0
  
  # Create X features
  X <- matrix(rnorm(N*5), nrow=N, ncol=5)
  
  # Parameters
  levels <- 1
  kNN <- 3
  embed_dim <- 4
  filter_power <- 2
  sigma <- 1e-3
  
  # Precompute A_feat from X using the same logic as in .fuse_graph()
  # We'll replicate a simplified version of that code here:
  all_knn <- hnsw_knn(X, k = kNN+1, distance = "l2")
  items <- all_knn$idx[, -1, drop=FALSE]
  
  norms <- sqrt(rowSums(X * X) + 1e-12)
  
  i_idx <- integer(N*kNN)
  j_idx <- integer(N*kNN)
  vals <- numeric(N*kNN)
  
  idx_ptr <- 1
  for (i in seq_len(N)) {
    nbrs <- items[i,]
    # It's possible some nodes have fewer neighbors if kNN is large relative to graph size
    nbrs <- nbrs[nbrs != 0]
    if (length(nbrs) == 0) next
    
    Xi <- X[i,,drop=FALSE]
    Xnbrs <- X[nbrs,,drop=FALSE]
    dotp <- Xi %*% t(Xnbrs)
    cos_sim <- as.numeric(dotp) / (norms[i]*norms[nbrs])
    
    len_nbrs <- length(nbrs)
    i_idx[idx_ptr:(idx_ptr+len_nbrs-1)] <- i
    j_idx[idx_ptr:(idx_ptr+len_nbrs-1)] <- nbrs
    vals[idx_ptr:(idx_ptr+len_nbrs-1)] <- cos_sim
    idx_ptr <- idx_ptr+len_nbrs
  }
  
  actual_size <- idx_ptr - 1
  i_idx <- i_idx[1:actual_size]
  j_idx <- j_idx[1:actual_size]
  vals <- vals[1:actual_size]
  
  A_feat <- sparseMatrix(i=i_idx, j=j_idx, x=vals, dims=c(N,N))
  A_feat <- (A_feat + t(A_feat))/2
  
  # Run graph_zoom with X (no A_feat), set seed for reproducibility
  set.seed(123)
  result_X <- graph_zoom(
    A = A, X = X, A_feat = NULL,
    levels = levels, kNN = kNN, embed_dim = embed_dim,
    filter_power = filter_power, sigma = sigma, verbose = FALSE
  )
  
  # Run graph_zoom with A_feat (no X), same seed
  set.seed(123)
  result_Afeat <- graph_zoom(
    A = A, X = NULL, A_feat = A_feat,
    levels = levels, kNN = kNN, embed_dim = embed_dim,
    filter_power = filter_power, sigma = sigma, verbose = FALSE
  )
  
  # Check that the final embeddings match
  # We'll use expect_equal with a tolerance in case of floating-point differences:
  expect_equal(result_X$E_original, result_Afeat$E_original, tolerance=1e-12)
  
  # Check that the coarsest embeddings match
  expect_equal(result_X$E_coarse, result_Afeat$E_coarse, tolerance=1e-12)
  
  # Check that the graph sequences match
  # Graph equality can be tricky due to floating point, but let's just check structure:
  expect_equal(result_X$graphs, result_Afeat$graphs, tolerance=1e-12)
  
  # Check that H_list matches
  expect_equal(result_X$H_list, result_Afeat$H_list, tolerance=1e-12)
})