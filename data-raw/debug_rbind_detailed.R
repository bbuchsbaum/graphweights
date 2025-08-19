# Debug the specific rbind issue when there are valid matches

library(Matrix)

# Recreate the exact scenario from threshold_adjacency
adj <- sparseMatrix(i = c(1, 1, 1, 2, 2, 3), j = c(1, 2, 3, 2, 3, 3), 
                   x = c(0.1, 0.5, 0.8, 0.3, 0.9, 0.2), dims = c(3, 3))

# What threshold_adjacency does first
k <- 2
rows <- parallel::mclapply(seq_len(nrow(adj)), mc.cores = 1, function(i) {
  ord <- order(adj[i, ], decreasing = TRUE)[seq_len(k)]
  cbind(i = i, j = ord, x = adj[i, ord])
})
rows <- do.call(rbind, rows)

m <- Matrix::sparseMatrix(i = rows[, 1], j = rows[, 2], x = rows[, 3], dims = dim(adj))

cat("Matrix after threshold_adjacency processing:\n")
print(m)

# Now debug psparse step by step
M <- as(m, "TsparseMatrix")
tri <- cbind(i = M@i + 1L, j = M@j + 1L, x = M@x)
tri_sym <- tri[tri[, 1] < tri[, 2], , drop = FALSE]

cat("\ntri:\n")
print(tri)
cat("tri_sym:\n")
print(tri_sym)

if (nrow(tri_sym) > 0) {
  tri_lookup <- paste(tri[, 1], tri[, 2], sep = ",")
  target_lookup <- paste(tri_sym[, 2], tri_sym[, 1], sep = ",")
  
  match_indices <- match(target_lookup, tri_lookup)
  valid_matches <- !is.na(match_indices)
  
  cat("match_indices:", match_indices, "\n")
  cat("valid_matches:", valid_matches, "\n")
  
  if (sum(valid_matches) > 0) {
    valid_tri_sym <- tri_sym[valid_matches, , drop = FALSE]
    valid_match_indices <- match_indices[valid_matches]
    
    cat("valid_tri_sym:\n")
    print(valid_tri_sym)
    cat("valid_match_indices:", valid_match_indices, "\n")
    
    new_x <- pmax(valid_tri_sym[, 3], M@x[valid_match_indices])
    
    out_tri <- rbind(
      cbind(valid_tri_sym[, 1:2], x = new_x),
      cbind(valid_tri_sym[, 2:1], x = new_x)
    )
    
    cat("out_tri:\n")
    print(out_tri)
    cat("out_tri structure:\n")
    str(out_tri)
    
    # The problematic part
    processed_positions <- c(which(tri[,1] < tri[,2])[valid_matches], valid_match_indices)
    unprocessed_mask <- !(1:nrow(tri) %in% processed_positions)
    diag_and_unmatched <- tri[unprocessed_mask, , drop = FALSE]
    
    cat("processed_positions:", processed_positions, "\n")
    cat("unprocessed_mask:", unprocessed_mask, "\n")
    cat("diag_and_unmatched:\n")
    print(diag_and_unmatched)
    cat("diag_and_unmatched structure:\n")
    str(diag_and_unmatched)
    
    # Check dimensions before rbind
    cat("out_tri dims:", dim(out_tri), "\n")
    cat("diag_and_unmatched dims:", dim(diag_and_unmatched), "\n")
    
    if (nrow(diag_and_unmatched) > 0) {
      cat("Checking column names:\n")
      cat("out_tri colnames:", colnames(out_tri), "\n")
      cat("diag_and_unmatched colnames:", colnames(diag_and_unmatched), "\n")
      
      # Try the rbind
      tryCatch({
        final_tri <- rbind(out_tri, diag_and_unmatched)
        cat("rbind succeeded!\n")
      }, error = function(e) {
        cat("rbind failed:", e$message, "\n")
        
        # Try to understand the structure difference
        cat("out_tri[1:2,]:\n")
        print(out_tri[1:2,])
        cat("diag_and_unmatched[1,]:\n")
        print(diag_and_unmatched[1,])
      })
    }
  }
}