# Test the fixed psparse function in isolation

library(Matrix)

# My fixed psparse function
psparse_fixed <- function(M, FUN, return_triplet = FALSE) {
  if (!inherits(M, "dgTMatrix")) M <- as(M, "TsparseMatrix")

  tri <- cbind(i = M@i + 1L, j = M@j + 1L, x = M@x)
  tri_sym <- tri[tri[, 1] < tri[, 2], , drop = FALSE]

  if (nrow(tri_sym) == 0) {
    # No off-diagonal elements to process
    if (return_triplet) return(tri)
    return(M)
  }

  # Fixed: Proper row-wise matching using vectorized lookup
  # Create lookup for (i,j) -> position in tri
  tri_lookup <- paste(tri[, 1], tri[, 2], sep = ",")
  target_lookup <- paste(tri_sym[, 2], tri_sym[, 1], sep = ",")  # Reversed pairs
  
  # Find positions of reversed pairs
  match_indices <- match(target_lookup, tri_lookup)
  
  # Handle case where some reversed pairs don't exist (asymmetric matrix)
  valid_matches <- !is.na(match_indices)
  
  cat("Debug: valid_matches =", valid_matches, "\n")
  cat("Debug: sum(valid_matches) =", sum(valid_matches), "\n")
  
  if (sum(valid_matches) == 0) {
    # No symmetric pairs found - return original matrix
    cat("Debug: Returning original matrix\n")
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
  processed_positions <- c(which(tri[,1] < tri[,2])[valid_matches], valid_match_indices)
  unprocessed_mask <- !(1:nrow(tri) %in% processed_positions)
  diag_and_unmatched <- tri[unprocessed_mask, , drop = FALSE]
  
  if (nrow(diag_and_unmatched) > 0) {
    out_tri <- rbind(out_tri, diag_and_unmatched)
  }

  if (return_triplet) return(out_tri)

  Matrix::sparseMatrix(i = out_tri[,1], j = out_tri[,2], x = out_tri[,3], dims = dim(M))
}

# Test with the problematic matrix
adj <- sparseMatrix(i = c(1, 1, 1, 2, 2, 3), j = c(1, 2, 3, 2, 3, 3), 
                   x = c(0.1, 0.5, 0.8, 0.3, 0.9, 0.2), dims = c(3, 3))

cat("Testing psparse_fixed:\n")
result <- psparse_fixed(adj, pmax)
cat("Success!\n")
print(result)

# Now test if it's called from threshold_adjacency correctly
cat("\nNow testing if the issue is in threshold_adjacency itself...\n")

# Let me check what threshold_adjacency actually calls
cat("Checking what threshold_adjacency does:\n")

# Simplified threshold_adjacency 
threshold_adjacency_test <- function(A, k = 5, type = c("normal", "mutual"), ncores = 1) {
  type <- match.arg(type)
  rows <- parallel::mclapply(seq_len(nrow(A)), mc.cores = ncores, function(i) {
    ord <- order(A[i, ], decreasing = TRUE)[seq_len(k)]
    cbind(i = i, j = ord, x = A[i, ord])
  })
  rows <- do.call(rbind, rows)

  m <- Matrix::sparseMatrix(i = rows[, 1], j = rows[, 2], x = rows[, 3],
                            dims = dim(A))

  cat("Before calling psparse_fixed with type:", type, "\n")
  psparse_fixed(m, if (type == "normal") pmax else pmin)
}

result2 <- threshold_adjacency_test(adj, k = 2)
cat("threshold_adjacency_test succeeded!\n")
print(result2)