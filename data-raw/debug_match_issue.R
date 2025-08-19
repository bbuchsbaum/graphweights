# Deeper investigation of the match issue in psparse function

library(Matrix)

# Create the same test case
adj <- sparseMatrix(i = c(1, 1, 1, 2, 2, 3), j = c(1, 2, 3, 2, 3, 3), 
                   x = c(0.1, 0.5, 0.8, 0.3, 0.9, 0.2), dims = c(3, 3))

M <- as(adj, "dgTMatrix")
tri <- cbind(i = M@i + 1L, j = M@j + 1L, x = M@x)
tri_sym <- tri[tri[, 1] < tri[, 2], , drop = FALSE]

cat("tri_sym (symmetric pairs only):\n")
print(tri_sym)

cat("\nTrying to find reverse pairs...\n")

# The problem: match() doesn't work properly with matrices for this use case
# Let's see what match is actually doing:
match_result <- match(tri_sym[, 2:1, drop = FALSE], tri[, 1:2, drop = FALSE])
cat("Match result:", match_result, "\n")

# This gives: 3 6 6 1 1 3, which is wrong!
# We need positions where (j,i) matches (i,j) in the original tri

cat("\nCorrect approach - manual matching:\n")
correct_indices <- integer(nrow(tri_sym))

for(k in 1:nrow(tri_sym)) {
  # Look for row where (j,i) = (tri_sym[k,2], tri_sym[k,1])
  target_i <- tri_sym[k, 2]  # j becomes i
  target_j <- tri_sym[k, 1]  # i becomes j
  
  for(m in 1:nrow(tri)) {
    if(tri[m, 1] == target_i && tri[m, 2] == target_j) {
      correct_indices[k] <- m
      cat("tri_sym[", k, "] = (", tri_sym[k,1], ",", tri_sym[k,2], 
          ") -> looking for (", target_i, ",", target_j, ") -> found at position", m, "\n")
      break
    }
  }
}

cat("Correct indices:", correct_indices, "\n")

# But some of these pairs don't exist in the original matrix!
# (2,1), (3,1), (3,2) are not in the original triplet list
# This means the matrix is not symmetric

cat("\nOriginal matrix symmetry check:\n")
print("adj == t(adj):")
print(as.matrix(adj == t(adj)))

cat("\nThe issue: The matrix is NOT symmetric!\n")
cat("But psparse assumes it's working with symmetric elements.\n")

# Let's check what happens when we use the actual values that exist
cat("\nValues that would be used with wrong match:\n")
wrong_values <- M@x[match_result]
cat("M@x[match_result]:", wrong_values, "\n")

# The problem is that match() with matrices is vectorizing incorrectly
# Let's see what happens step by step

cat("\nDebug match() behavior:\n")
target_matrix <- tri_sym[, 2:1, drop = FALSE]
search_matrix <- tri[, 1:2, drop = FALSE]

cat("target_matrix:\n")
print(target_matrix)
cat("search_matrix:\n") 
print(search_matrix)

# match() treats these as vectors, not row-wise matching!
cat("as.vector(target_matrix):", as.vector(target_matrix), "\n")
cat("as.vector(search_matrix):", as.vector(search_matrix), "\n")

match_vector_result <- match(as.vector(target_matrix), as.vector(search_matrix))
cat("match(as.vector(target), as.vector(search)):", match_vector_result, "\n")

# This explains the dimension mismatch warning!
# The result has length 6 (2 cols * 3 rows) but we're trying to use it 
# to index M@x which has length 6, but in a different context