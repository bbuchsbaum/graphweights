# Debug threshold_adjacency dimension mismatch issue

library(Matrix)
source("../R/knn_weights.R")

# Reproduce the exact test case
adj <- sparseMatrix(i = c(1, 1, 1, 2, 2, 3), j = c(1, 2, 3, 2, 3, 3), 
                   x = c(0.1, 0.5, 0.8, 0.3, 0.9, 0.2), dims = c(3, 3))

cat("Original adjacency matrix:\n")
print(adj)

cat("\nMatrix structure:\n")
print(str(adj))

cat("\nConverting to dgTMatrix:\n")
M <- as(adj, "dgTMatrix")
print(M)

cat("\nExtracting triplet format:\n")
tri <- cbind(i = M@i + 1L, j = M@j + 1L, x = M@x)
print(tri)

cat("\nFiltering symmetric part (i < j):\n")
tri_sym <- tri[tri[, 1] < tri[, 2], , drop = FALSE]
print(tri_sym)

cat("\nDimensions:\n")
cat("tri_sym dimensions:", dim(tri_sym), "\n")
cat("tri_sym[, 2:1] dimensions:", dim(tri_sym[, 2:1, drop = FALSE]), "\n")
cat("tri[, 1:2] dimensions:", dim(tri[, 1:2, drop = FALSE]), "\n")

# The problematic match operation
cat("\nTrying the match operation:\n")
cat("tri_sym[, 2:1]:\n")
print(tri_sym[, 2:1, drop = FALSE])

cat("\ntri[, 1:2]:\n")  
print(tri[, 1:2, drop = FALSE])

# This is where the error occurs
tryCatch({
  match_result <- match(tri_sym[, 2:1, drop = FALSE], tri[, 1:2, drop = FALSE])
  cat("Match result:", match_result, "\n")
}, error = function(e) {
  cat("Error in match:", e$message, "\n")
})

# Let's see what M@x looks like
cat("\nM@x values:", M@x, "\n")

# Try to understand the issue by manually debugging
cat("\nDebugging the match issue:\n")

# The match is trying to find positions of reversed pairs
for(i in 1:nrow(tri_sym)) {
  reversed_pair <- tri_sym[i, 2:1]
  cat("Looking for reversed pair [", reversed_pair[1], ",", reversed_pair[2], "] in tri\n")
  
  # Find manually
  for(j in 1:nrow(tri)) {
    if(tri[j, 1] == reversed_pair[1] && tri[j, 2] == reversed_pair[2]) {
      cat("  Found at position", j, "with value", tri[j, 3], "\n")
    }
  }
}

# The real test - call the actual function
cat("\nCalling threshold_adjacency:\n")
tryCatch({
  result <- threshold_adjacency(adj, k = 2)
  cat("Success! Result:\n")
  print(result)
}, error = function(e) {
  cat("Error:", e$message, "\n")
}, warning = function(w) {
  cat("Warning:", w$message, "\n")
})