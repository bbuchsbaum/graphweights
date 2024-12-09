#' @title Multi-level Graph Coarsening by Local Variation (Edge-based)
#' @description Implements the multi-level edge-based local variation coarsening algorithm 
#' from Loukas (2019), "Graph Reduction with Spectral and Cut Guarantees" (JMLR).
#'
#' Given a sparse adjacency matrix of an undirected, weighted graph and a target dimension,
#' this algorithm produces a reduced graph that approximates a given spectral subspace 
#' (e.g., the first k eigenvectors of the Laplacian) in the restricted spectral approximation sense.
#'
#' The algorithm follows Section 4 of the referenced paper and uses the local variation approach 
#' to ensure good spectral approximation. It performs multiple levels of coarsening, at each level 
#' selecting candidate sets (edges) that minimize local variation costs (Sections 4.1-4.3).
#'
#' It also leverages Proposition 17 to update the error estimate and decides when to stop based 
#' on the given threshold epsilon_prime.
#'
#' @param A A sparse adjacency matrix of class "dgCMatrix" representing an undirected graph.
#' @param n Target number of vertices after coarsening.
#' @param k Dimension of the subspace to preserve (e.g., first k eigenvectors).
#' @param epsilon_prime Error threshold for restricted spectral approximation.
#' @param max_levels Maximum number of coarsening levels (default: \code{ceiling(log2(nrow(A)/n))}).
#' @param verbose If TRUE, prints progress messages.
#'
#' @return A list with elements:
#' \item{A_coarse}{The adjacency matrix of the coarsened graph.}
#' \item{P}{The overall reduction matrix \(P\).}
#' \item{info}{A list with details about intermediate steps, final approximation error, and levels used.}
#'
#' @references 
#' Loukas, A. (2019). "Graph Reduction with Spectral and Cut Guarantees." 
#' Journal of Machine Learning Research, 20(116):1–42.
#' 
#' @export
coarsen_graph <- function(A, n, k=10, epsilon_prime=0.5, max_levels=NULL, verbose=FALSE) {
  if(!inherits(A, "dgCMatrix")) stop("A must be a sparse 'dgCMatrix'.")
  N <- nrow(A)
  if(is.null(max_levels)) {
    max_levels <- ceiling(log2(N/n))
  }
  if(k >= N) stop("k must be smaller than the number of vertices.")
  
  # Compute Laplacian L = D - A
  d <- Matrix::rowSums(A)
  L <- Matrix::Diagonal(x=d) - A
  
  # Compute first k eigenvectors of L for the chosen subspace R = span(U_k)
  # (Section 3.2 of paper)
  if(verbose) message("Computing first ", k, " eigenvectors of the Laplacian...")
  eig <- RSpectra::eigs(L, k = k, which = "SM")
  U_k <- eig$vectors
  lambda <- eig$values
  
  # Handle zero or near-zero eigenvalues by adding a small offset for stability
  offset <- 1e-12
  lambda_reg <- lambda + offset
  inv_sqrt_lambda <- 1 / sqrt(lambda_reg) 
  # B_0 = U_k * Lambda_k^{+1/2} (Remark in Section 4.1)
  B_prev <- U_k %*% diag(inv_sqrt_lambda, k, k)
  
  # Initialize multi-level coarsening
  L_cur <- L
  N_cur <- N
  epsilon_cur <- 0
  
  # We store transformations:
  # P_all: list of P at each level, start with identity
  P_all <- list(Matrix::Diagonal(N))
  
  level <- 0
  while(N_cur > n && epsilon_cur < epsilon_prime && level < max_levels) {
    level <- level + 1
    if(verbose) message("Level ", level, ": current size=", N_cur, ", epsilon=", epsilon_cur)
    
    # sigma_prime for this level (Proposition 17)
    sigma_prime <- ((1+epsilon_prime)/(1+epsilon_cur))-1
    
    # Single-level coarsening (Algorithm 2)
    # We provide B_prev, the subspace from previous level
    # According to Section 4.1, after coarsening, we must update B for next level.
    # B_{ell} = P_ell * B_{ell-1}, and then we must re-orthonormalize and 
    # apply (B_{ell}^T L_{ell} B_{ell})^{+1/2} inverse if needed.
    # However, the local variation cost uses B_{ell-1}. We'll compute after finalizing coarsening.
    
    res_level <- coarsen_one_level(L_cur, B_prev, n = max(n, floor(N_cur/2)), sigma_prime = sigma_prime, verbose = verbose)
    L_new <- res_level$L_new
    P_level <- res_level$P
    sigma_level <- res_level$sigma
    
    # Update epsilon (Proposition 17)
    epsilon_cur <- (1+epsilon_cur)*(1+sigma_level)-1
    
    # Update L and N
    L_cur <- L_new
    N_cur <- nrow(L_new)
    
    # Update overall P: P = P * P_level
    # We only store P_level and will multiply at the end for final P,
    # but to keep track, we do it step by step:
    # Instead of recomputing a huge product each time, we store cumulative transformations.
    P_all[[length(P_all)+1]] <- P_level
    
    # Update B for next iteration:
    # B_ell = P_ell * B_{ell-1}, from Section 4.1:
    B_curr <- P_level %*% B_prev
    
    # Now we must apply (B_curr^T L_cur B_curr)^{+1/2} pseudo-inverse to maintain consistency.
    # This ensures that at each level we have B_ell = P_ell ... P_1 U_k Lambda_k^{+1/2}
    # This keeps B_ell suitable for variation cost computations on subsequent levels.
    M_mat <- t(B_curr) %*% L_cur %*% B_curr
    # Compute pseudo-inverse of M_mat^{1/2}
    # eigen-decompose M_mat:
    M_eig <- eigen(M_mat, symmetric=TRUE)
    # Handle small eigenvalues:
    M_eig_val <- M_eig$values
    M_eig_val_reg <- M_eig_val + offset
    inv_sqrt_M <- M_eig$vectors %*% diag(1/sqrt(M_eig_val_reg), k, k) %*% t(M_eig$vectors)
    # B_ell = B_curr * inv_sqrt_M
    B_prev <- B_curr %*% inv_sqrt_M
    
    if(verbose) message("After level ", level, ": size=", N_cur, ", epsilon=", epsilon_cur)
  }
  
  # Combine all P's to get final P:
  # P = P_c ... P_1
  # We stored P_all: P_all[[1]] = I, P_all[[2]]=P_1, ..., P_all[[end]]=P_c
  # Compute final P carefully to avoid numerical instability.
  if(verbose) message("Combining all coarsening matrices P...")
  P_final <- P_all[[1]]
  for(i in 2:length(P_all)) {
    P_final <- P_final %*% P_all[[i]]
  }
  
  # Construct A_coarse from L_cur:
  d_c <- Matrix::diag(L_cur)
  A_coarse <- Matrix::Diagonal(x=d_c) - L_cur
  
  # Return results
  list(
    A_coarse = A_coarse,
    P = P_final,
    info = list(
      final_epsilon = epsilon_cur,
      final_n = N_cur,
      levels = level
    )
  )
}


#' @title Single-level Coarsening by Local Variation (Edge-based)
#' @description Implements one step of the local variation coarsening (Algorithm 2 in Loukas (2019)).
#'
#' Given a Laplacian L and subspace B from the previous level, a target size n, and a variation threshold sigma_prime,
#' this function returns a coarsened Laplacian L_new and the coarsening matrix P that achieve 
#' a variation cost sigma <= sigma_prime if possible.
#'
#' @param L The combinatorial Laplacian at the current level.
#' @param B The subspace matrix B_{ell-1}, as described in Section 4.1 of the paper.
#' @param n Target size after coarsening.
#' @param sigma_prime Variation threshold for this level.
#' @param verbose If TRUE, prints progress messages.
#'
#' @return A list:
#' \item{L_new}{The coarsened Laplacian.}
#' \item{P}{The coarsening matrix at this level.}
#' \item{sigma}{The variation cost achieved at this level.}
#'
#' @references 
#' Loukas, A. (2019) "Graph Reduction with Spectral and Cut Guarantees." JMLR.
#'
#' @export
coarsen_one_level <- function(L, B, n, sigma_prime, verbose=FALSE) {
  N <- nrow(L)
  if(n >= N) {
    # No reduction needed
    return(list(L_new=L, P=Matrix::Diagonal(N), sigma=0))
  }
  
  # From L, get A and degrees for edge-based candidates
  d <- Matrix::diag(L)
  A <- Matrix::Diagonal(x=d) - L
  
  # Extract edges i<j from A
  coo <- summary(A)
  edges <- coo[coo$i < coo$j & coo$x > 0,]
  # edges now hold i,j, and x
  # We'll consider candidate sets C={i,j}.
  # Variation cost for set C=(i,j) is defined by Equation (6) in the paper:
  # cost(C) = || Pi_C^{\perp} B ||_{L_C}^2 / (|C|-1)
  # |C|=2, so denominator=1.
  #
  # Pi_C^{\perp} projects onto the subspace orthogonal to the constant vector on C.
  # For C={i,j}, Pi_C^{\perp} applied to any vector b gives [ (b(i)-b(j))/2; -(b(i)-b(j))/2 ].
  #
  # L_C is defined in Section 4.2 (Equation (5)):
  # For edges:
  # W_C(i,j)=W(i,j) inside C;
  # Edges going outside C are counted twice.
  
  # Precompute row sums once for efficiency
  row_sums_A <- Matrix::rowSums(A)
  
  compute_edge_cost <- function(i, j, B, L, A, row_sums_A) {
    # Outside edges from i: sum(A[i,]) - A[i,i] - A[i,j]
    w_ij <- A[i,j]
    out_i <- row_sums_A[i] - A[i,i] - w_ij
    out_j <- row_sums_A[j] - A[j,j] - w_ij
    
    deg_iC <- w_ij + 2*out_i
    deg_jC <- w_ij + 2*out_j
    # L_C:
    # [deg_iC, -w_ij; -w_ij, deg_jC]
    L_C <- matrix(c(deg_iC, -w_ij, -w_ij, deg_jC), 2, 2)
    
    # Pi_C^perp B: For each column b of B:
    # delta_l = (b(i)-b(j))/2
    # vector = [delta_l; -delta_l]
    b_i <- B[i,]
    b_j <- B[j,]
    delta_vec <- (b_i - b_j)/2
    
    # cost = sum over l: v_l^T L_C v_l with v_l=[delta_l; -delta_l]
    # v_l = delta_l * [1; -1]
    # Let's do a vectorized compute:
    # cost_col_l = (delta_l^2) * [1, -1] * L_C * [1; -1]
    # Compute q = L_C*[1; -1]
    q <- L_C %*% c(1,-1)
    # q = [deg_iC + w_ij; -w_ij - deg_jC] actually
    # cost per column l:
    # cost_l = delta_l^2 * ( [1,-1] %*% q )
    # [1,-1]%*%q = q[1]*1 + q[2]*(-1) = deg_iC+w_ij -(-w_ij - deg_jC)
    # = deg_iC + w_ij + w_ij + deg_jC
    # = deg_iC + deg_jC + 2*w_ij
    # Wait, check correctness:
    # Actually, v_l^T L_C v_l = delta_l^2 ([1,-1]L_C[1;-1])
    # = delta_l^2 ([1,-1]*q)
    # q = L_C*[1;-1]
    # Let's just compute it directly once:
    val <- c(1,-1) %*% q
    # val = deg_iC + deg_jC + 2*w_ij
    
    cost_cols <- (delta_vec^2)*val
    cost_val <- sum(cost_cols)
    return(cost_val)
  }
  
  if(verbose) message("Computing initial costs for edges...")
  # Precompute all edge costs:
  edge_costs <- apply(edges, 1, function(e) {
    i <- e["i"]; j <- e["j"]
    compute_edge_cost(i,j,B,L,A,row_sums_A)
  })
  
  ord <- order(edge_costs)
  edges_sorted <- cbind(edges[ord,], cost=edge_costs[ord])
  
  # Algorithm 2 loop
  marked <- rep(FALSE, N)
  sigma_ell_sq <- 0
  N_ell <- N
  Psets <- list()
  
  idx <- 1
  # Iterate over candidates in ascending cost order
  while(idx <= nrow(edges_sorted) && N_ell > n && sqrt(sigma_ell_sq) <= sigma_prime) {
    e <- edges_sorted[idx,]
    i <- e["i"]; j <- e["j"]
    s <- e["cost"]
    if(!marked[i] && !marked[j]) {
      candidate_sigma_sq <- sigma_ell_sq + s
      if(sqrt(candidate_sigma_sq) <= sigma_prime) {
        # Accept this contraction
        marked[i] <- TRUE
        marked[j] <- TRUE
        Psets[[length(Psets)+1]] <- c(i,j)
        sigma_ell_sq <- candidate_sigma_sq
        N_ell <- N_ell - 1
      }
    }
    idx <- idx + 1
  }
  
  # Add singletons for unmarked vertices
  singletons <- which(!marked)
  
  # Construct P_ell
  total_new_vertices <- length(Psets) + length(singletons)
  P <- Matrix::Matrix(0, nrow=total_new_vertices, ncol=N, sparse=TRUE)
  
  # Contracted sets (edges)
  for(r in seq_along(Psets)) {
    C <- Psets[[r]]
    # uniform weights = 1/2 for edges
    P[r,C] <- 1/2
  }
  
  # Singletons
  for(r in seq_along(singletons)) {
    P[length(Psets)+r, singletons[r]] <- 1
  }
  
  # Compute L_new = P^{\mp} L P^{+} (Proposition 6,7)
  # For coarsening:
  # P^{+} = P^{T} D_P^{-2}, D_P(r,r)=||P(r,:)||_2
  row_norms <- sqrt(Matrix::rowSums(P^2))
  invD2 <- Matrix::Diagonal(x=1/(row_norms^2))
  
  P_plus <- t(P) %*% invD2
  P_mp <- invD2 %*% P
  L_new <- P_mp %*% L %*% P_plus
  
  sigma_ell <- sqrt(sigma_ell_sq)
  
  if(verbose) {
    message("Single-level coarsening done. Achieved size: ", nrow(L_new),
            " sigma=", sigma_ell)
  }
  
  list(L_new=L_new, P=P, sigma=sigma_ell)
}