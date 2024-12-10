#' @title Multi-level Graph Coarsening by Local Variation (Edge-based)
#' @description Implements the multi-level edge-based local variation coarsening 
#' based on Loukas (2019). The method updates the Laplacian and a subspace U at each level,
#' ensuring stable re-orthonormalization and forcing at least one contraction if needed.
#'
#' @param A A sparse adjacency matrix (dgCMatrix) of the undirected graph.
#' @param n Target number of vertices after coarsening (must be < N).
#' @param k Dimension of the subspace to preserve (e.g., first k eigenvectors of L).
#' @param epsilon_prime Error threshold for restricted spectral approximation (used as a heuristic).
#' @param max_levels Maximum levels allowed (default: log2(N/n)), just to prevent infinite loops.
#' @param verbose If TRUE, prints progress messages.
#'
#' @return A list:
#' \item{A_coarse}{The coarsened adjacency matrix.}
#' \item{P}{The overall reduction matrix mapping original N to final n.}
#' \item{info}{A list with final_epsilon, final_n, levels.}
#'
#' @references Loukas, A. (2019). "Graph Reduction with Spectral and Cut Guarantees."
#' Journal of Machine Learning Research, 20(116):1–42.
#' 
#' @export
#' @importFrom RSpectra eigs
#' @importFrom Matrix Diagonal 
coarsen_graph_multilevel <- function(A, n, k=2, epsilon_prime=0.5, max_levels=NULL, verbose=FALSE) {
  if(!inherits(A, "dgCMatrix")) stop("A must be a 'dgCMatrix'.")
  N <- nrow(A)
  if(n >= N) stop("n must be less than N.")

  if(is.null(max_levels)) max_levels <- ceiling(log2(N/n)) + 5

  # Compute initial Laplacian and subspace
  d <- rowSums(A)
  L <- Diagonal(x=d) - A
  # Compute first k eigenvectors of L
  if(verbose) message("Computing first ", k, " eigenvectors of the Laplacian...")
  eig <- RSpectra::eigs(L, k=k, which="SM")
  U <- eig$vectors
  lambda <- eig$values

  offset <- 1e-9
  L_cur <- L
  N_cur <- N
  epsilon_cur <- 0
  level <- 0
  P_all <- list()

  while(N_cur > n && level < max_levels) {
    level <- level + 1
    if(verbose) message("Level ", level, ": current size=", N_cur, ", epsilon=", epsilon_cur)

    # sigma_prime heuristic:
    sigma_prime <- epsilon_prime + 1e-12

    # Perform one level of coarsening
    res_level <- coarsen_one_level_localvar(L_cur, U, n, sigma_prime, force_contraction=FALSE, verbose=verbose)
    if(res_level$N_new == N_cur && N_cur > n) {
      # No reduction, force contraction
      if(verbose) message("No reduction achieved. Forcing contraction.")
      res_level <- coarsen_one_level_localvar(L_cur, U, n, sigma_prime, force_contraction=TRUE, verbose=verbose)
    }

    if(res_level$N_new == N_cur) {
      # Still no reduction possible
      if(verbose) message("Unable to reduce further, stopping.")
      break
    }

    L_cur <- res_level$L_new
    U <- res_level$U_new  # Updated U
    N_cur <- res_level$N_new
    P_ell <- res_level$P
    sigma_level <- res_level$sigma
    epsilon_cur <- epsilon_cur + sigma_level

    P_all[[level]] <- P_ell

    # Re-orthonormalize U to ensure stability
    M_mat <- t(U) %*% L_cur %*% U
    M_eig <- eigen(M_mat, symmetric=TRUE)
    M_val <- M_eig$values
    min_val <- min(M_val)
    if(min_val <= 0) {
      shift <- -min_val + offset
      M_val_reg <- M_val + shift
    } else {
      M_val_reg <- M_val + offset
    }

    inv_sqrt_M <- M_eig$vectors %*% diag(1/sqrt(M_val_reg), k, k) %*% t(M_eig$vectors)
    U <- U %*% inv_sqrt_M

    if(verbose) message("After level ", level, ": size=", N_cur, ", epsilon=", epsilon_cur)
  }

  # Combine all P to get final P
  if(verbose) message("Combining all coarsening matrices P...")
  if(length(P_all) == 0) {
    # No coarsening done
    P_final <- Diagonal(N)
    A_coarse <- A
  } else {
    # P_all[[1]] maps R^N -> R^{N_1}
    # P_all[[2]] maps R^{N_1} -> R^{N_2}, ...
    # Final P = P_c * ... * P_1
    P_final <- P_all[[length(P_all)]]
    if(length(P_all) > 1) {
      for(i in (length(P_all)-1):1) {
        # Check dimensions before multiplication
        if(ncol(P_final) != nrow(P_all[[i]])) {
          stop("Dimension mismatch when combining P matrices.")
        }
        P_final <- P_final %*% P_all[[i]]
      }
    }

    d_c <- diag(L_cur)
    A_coarse <- Diagonal(x=d_c) - L_cur
  }

  list(
    A_coarse = A_coarse,
    P = P_final,
    info = list(
      final_epsilon = epsilon_cur,
      final_n = nrow(A_coarse),
      levels = level
    )
  )
}


#' @title Single-level local variation coarsening step (Edge-based)
#' @description Performs one coarsening step:
#' - Computes edge costs
#' - Selects edges to contract (if none meet threshold, force one contraction if force_contraction=TRUE)
#' - Constructs P, updates L and U
#' 
#' @param L Current Laplacian.
#' @param U Current subspace U (N x k).
#' @param n Target size after coarsening (N_new <= n ideally).
#' @param sigma_prime Variation threshold.
#' @param force_contraction If TRUE, force at least one edge contraction if none selected normally.
#' @param verbose If TRUE, prints messages.
#'
#' @return A list:
#' \item{L_new}{Reduced Laplacian.}
#' \item{P}{Coarsening matrix.}
#' \item{U_new}{Updated subspace.}
#' \item{N_new}{New dimension.}
#' \item{sigma}{Approx. variation cost (for epsilon update).}
#'
#' @export
coarsen_one_level_localvar <- function(L, U, n, sigma_prime, force_contraction=FALSE, verbose=FALSE) {
  N <- nrow(L)
  if(n >= N) {
    return(list(L_new=L, P=Diagonal(N), U_new=U, N_new=N, sigma=0))
  }

  d <- diag(L)
  A <- Diagonal(x=d) - L

  coo <- summary(A)
  edges <- coo[coo$i < coo$j & coo$x > 0,]
  if(nrow(edges)==0) {
    # No edges, cannot reduce
    if(verbose) message("No edges to contract.")
    return(list(L_new=L, P=Diagonal(N), U_new=U, N_new=N, sigma=0))
  }

  row_sums_A <- rowSums(A)
  compute_edge_cost <- function(i, j, U, L, A, row_sums_A) {
    # Variation cost as in previous attempts
    w_ij <- A[i,j]
    out_i <- row_sums_A[i] - A[i,i] - w_ij
    out_j <- row_sums_A[j] - A[j,j] - w_ij
    deg_iC <- w_ij + 2*out_i
    deg_jC <- w_ij + 2*out_j
    L_C <- matrix(c(deg_iC, -w_ij, -w_ij, deg_jC),2,2)

    b_i <- U[i,]
    b_j <- U[j,]
    delta_vec <- (b_i - b_j)/2
    q <- L_C %*% c(1,-1)
    val <- as.numeric(c(1,-1) %*% q)
    cost_cols <- (delta_vec^2)*val
    sum(cost_cols)
  }

  # Compute costs
  costs <- apply(edges, 1, function(e) {
    i <- as.integer(e["i"])
    j <- as.integer(e["j"])
    compute_edge_cost(i,j,U,L,A,row_sums_A)
  })
  edges <- cbind(edges, cost=costs)

  # We want to contract edges that keep sqrt(sigma_ell_sq) <= sigma_prime
  # We'll try greedily: sort by cost, pick edges that keep the condition
  edges <- edges[order(edges[,"cost"]),]

  marked <- rep(FALSE, N)
  sigma_ell_sq <- 0
  N_ell <- N
  Psets <- list()

  # Attempt normal contraction with sigma_prime
  idx <- 1
  while(idx <= nrow(edges) && N_ell > n) {
    e <- edges[idx,]
    i <- as.integer(e["i"])
    j <- as.integer(e["j"])
    s <- as.numeric(e["cost"])
    if(!marked[i] && !marked[j]) {
      candidate_sigma_sq <- sigma_ell_sq + s
      if(sqrt(candidate_sigma_sq) <= sigma_prime) {
        # Accept contraction
        marked[i] <- TRUE
        marked[j] <- TRUE
        Psets[[length(Psets)+1]] <- c(i,j)
        sigma_ell_sq <- candidate_sigma_sq
        N_ell <- N_ell - 1
      }
    }
    idx <- idx + 1
  }

  # If no contraction and force_contraction = TRUE
  if(N_ell == N && n < N && force_contraction) {
    # Force contract the cheapest edge not marked
    for(idx in seq_len(nrow(edges))) {
      e <- edges[idx,]
      i <- as.integer(e["i"])
      j <- as.integer(e["j"])
      s <- as.numeric(e["cost"])
      if(!marked[i] && !marked[j]) {
        marked[i] <- TRUE
        marked[j] <- TRUE
        Psets[[length(Psets)+1]] <- c(i,j)
        sigma_ell_sq <- sigma_ell_sq + s
        N_ell <- N_ell - 1
        break
      }
    }
  }

  if(N_ell == N) {
    # No edges contracted
    if(verbose) message("No edges contracted, returning unchanged.")
    return(list(L_new=L, P=Diagonal(N), U_new=U, N_new=N, sigma=0))
  }

  singletons <- which(!marked)
  total_new_vertices <- length(Psets) + length(singletons)

  P <- Matrix(0, nrow=total_new_vertices, ncol=N, sparse=TRUE)
  row_idx <- 1
  for(C in Psets) {
    P[row_idx, C] <- 1/2
    row_idx <- row_idx+1
  }
  for(vx in singletons) {
    P[row_idx, vx] <- 1
    row_idx <- row_idx+1
  }

  # Compute L_new
  row_norms <- sqrt(rowSums(P^2))
  invD2 <- Diagonal(x=1/(row_norms^2))
  P_plus <- t(P) %*% invD2
  P_mp <- invD2 %*% P
  L_new <- P_mp %*% L %*% P_plus

  # Update U: U_new = P U
  U_new <- P %*% U

  sigma_ell <- sqrt(sigma_ell_sq)

  if(verbose) {
    message("Single-level coarsening done. Achieved size: ", nrow(L_new),
            " sigma=", sigma_ell)
  }

  list(L_new=L_new, P=P, U_new=U_new, N_new=nrow(L_new), sigma=sigma_ell)
}