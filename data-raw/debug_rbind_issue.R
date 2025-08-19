# Debug the rbind column mismatch issue

library(Matrix)

# Test the problematic case step by step
adj <- sparseMatrix(i = c(1, 1, 1, 2, 2, 3), j = c(1, 2, 3, 2, 3, 3), 
                   x = c(0.1, 0.5, 0.8, 0.3, 0.9, 0.2), dims = c(3, 3))

M <- as(adj, "TsparseMatrix")
tri <- cbind(i = M@i + 1L, j = M@j + 1L, x = M@x)
tri_sym <- tri[tri[, 1] < tri[, 2], , drop = FALSE]

cat("tri:\n")
print(tri)
cat("Dimensions:", dim(tri), "\n\n")

cat("tri_sym:\n") 
print(tri_sym)
cat("Dimensions:", dim(tri_sym), "\n\n")

if (nrow(tri_sym) > 0) {
  tri_lookup <- paste(tri[, 1], tri[, 2], sep = ",")
  target_lookup <- paste(tri_sym[, 2], tri_sym[, 1], sep = ",")
  
  cat("tri_lookup:", tri_lookup, "\n")
  cat("target_lookup:", target_lookup, "\n")
  
  match_indices <- match(target_lookup, tri_lookup)
  cat("match_indices:", match_indices, "\n")
  
  valid_matches <- !is.na(match_indices)
  cat("valid_matches:", valid_matches, "\n")
  
  if (sum(valid_matches) > 0) {
    valid_tri_sym <- tri_sym[valid_matches, , drop = FALSE]
    valid_match_indices <- match_indices[valid_matches]
    
    cat("valid_tri_sym:\n")
    print(valid_tri_sym)
    cat("valid_match_indices:", valid_match_indices, "\n")
    
    # Try to simulate the FUN operation
    new_x <- pmax(valid_tri_sym[, 3], M@x[valid_match_indices])
    cat("new_x:", new_x, "\n")
    
    out_tri <- rbind(
      cbind(valid_tri_sym[, 1:2], x = new_x),
      cbind(valid_tri_sym[, 2:1], x = new_x)
    )
    
    cat("out_tri:\n")
    print(out_tri)
    cat("out_tri dimensions:", dim(out_tri), "\n")
    
    # Now the problematic part
    processed_positions <- c(which(tri[,1] < tri[,2])[valid_matches], valid_match_indices)
    cat("processed_positions:", processed_positions, "\n")
    
    unprocessed_mask <- !(1:nrow(tri) %in% processed_positions)
    cat("unprocessed_mask:", unprocessed_mask, "\n")
    
    diag_and_unmatched <- tri[unprocessed_mask, , drop = FALSE]
    cat("diag_and_unmatched:\n")
    print(diag_and_unmatched)
    cat("diag_and_unmatched dimensions:", dim(diag_and_unmatched), "\n")
    
    if (nrow(diag_and_unmatched) > 0) {
      cat("Trying rbind...\n")
      cat("out_tri columns:", ncol(out_tri), "\n")
      cat("diag_and_unmatched columns:", ncol(diag_and_unmatched), "\n")
      
      if (ncol(out_tri) != ncol(diag_and_unmatched)) {
        cat("ERROR: Column mismatch!\n")
        cat("out_tri colnames:", colnames(out_tri), "\n")
        cat("diag_and_unmatched colnames:", colnames(diag_and_unmatched), "\n")
      }
    }
  }
}