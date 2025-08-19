# Fix for the psparse function

library(Matrix)

# Current broken implementation (simplified)
psparse_broken <- function(M, FUN, return_triplet = FALSE) {
  if (!inherits(M, "dgTMatrix")) M <- as(M, "dgTMatrix")
  
  tri <- cbind(i = M@i + 1L, j = M@j + 1L, x = M@x)
  tri_sym <- tri[tri[, 1] < tri[, 2], , drop = FALSE]
  
  # BROKEN: This match doesn't work correctly for row-wise matching
  new_x <- FUN(
    tri_sym[, 3],
    M@x[match(tri_sym[, 2:1, drop = FALSE], tri[, 1:2, drop = FALSE])]
  )
  
  out_tri <- rbind(
    cbind(tri_sym[, 1:2], x = new_x),
    cbind(tri_sym[, 2:1], x = new_x)
  )
  
  if (return_triplet) return(out_tri)
  
  # Assuming triplet_to_matrix exists
  # triplet_to_matrix(out_tri, dim(M))
  Matrix::sparseMatrix(i = out_tri[,1], j = out_tri[,2], x = out_tri[,3], dims = dim(M))
}

# Fixed implementation
psparse_fixed <- function(M, FUN, return_triplet = FALSE) {
  if (!inherits(M, "dgTMatrix")) M <- as(M, "TsparseMatrix")  # Use new syntax
  
  tri <- cbind(i = M@i + 1L, j = M@j + 1L, x = M@x)
  tri_sym <- tri[tri[, 1] < tri[, 2], , drop = FALSE]
  
  if (nrow(tri_sym) == 0) {
    # No off-diagonal elements to process
    if (return_triplet) return(tri)
    return(M)
  }
  
  # FIXED: Proper row-wise matching using vectorized approach
  # Create lookup for (i,j) -> position in tri
  tri_lookup <- paste(tri[, 1], tri[, 2], sep = ",")
  target_lookup <- paste(tri_sym[, 2], tri_sym[, 1], sep = ",")  # Reversed pairs
  
  # Find positions of reversed pairs
  match_indices <- match(target_lookup, tri_lookup)
  
  # Handle case where some reversed pairs don't exist (asymmetric matrix)
  valid_matches <- !is.na(match_indices)
  
  if (sum(valid_matches) == 0) {
    # No symmetric pairs found - return original matrix
    if (return_triplet) return(tri)
    return(M)
  }
  
  # Only process valid symmetric pairs
  valid_tri_sym <- tri_sym[valid_matches, , drop = FALSE]
  valid_match_indices <- match_indices[valid_matches]
  
  new_x <- FUN(
    valid_tri_sym[, 3],
    M@x[valid_match_indices]
  )
  
  out_tri <- rbind(
    cbind(valid_tri_sym[, 1:2], x = new_x),
    cbind(valid_tri_sym[, 2:1], x = new_x)
  )
  
  # Add back the diagonal and unmatched elements
  diag_and_unmatched <- tri[!(1:nrow(tri) %in% c(which(tri[,1] < tri[,2])[valid_matches], valid_match_indices)), , drop = FALSE]
  
  if (nrow(diag_and_unmatched) > 0) {
    out_tri <- rbind(out_tri, diag_and_unmatched)
  }
  
  if (return_triplet) return(out_tri)
  
  Matrix::sparseMatrix(i = out_tri[,1], j = out_tri[,2], x = out_tri[,3], dims = dim(M))
}

# Test with the problematic case
cat("=== Testing with problematic matrix ===\n")
adj <- sparseMatrix(i = c(1, 1, 1, 2, 2, 3), j = c(1, 2, 3, 2, 3, 3), 
                   x = c(0.1, 0.5, 0.8, 0.3, 0.9, 0.2), dims = c(3, 3))

cat("Original matrix:\n")
print(adj)

cat("\nTesting broken version:\n")
tryCatch({
  result_broken <- psparse_broken(adj, pmax)
  cat("Success (unexpected!):\n")
  print(result_broken)
}, warning = function(w) {
  cat("Warning:", w$message, "\n")
}, error = function(e) {
  cat("Error:", e$message, "\n")
})

cat("\nTesting fixed version:\n")
tryCatch({
  result_fixed <- psparse_fixed(adj, pmax)
  cat("Success:\n")
  print(result_fixed)
}, warning = function(w) {
  cat("Warning:", w$message, "\n")
}, error = function(e) {
  cat("Error:", e$message, "\n")
})

# Test with a symmetric matrix to make sure it still works
cat("\n=== Testing with symmetric matrix ===\n")
adj_sym <- adj + t(adj)
diag(adj_sym) <- diag(adj)  # Keep original diagonal
cat("Symmetric matrix:\n")
print(adj_sym)

cat("\nFixed version on symmetric matrix:\n")
result_sym <- psparse_fixed(adj_sym, pmax)
print(result_sym)