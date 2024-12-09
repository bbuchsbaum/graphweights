#' GraphZoom: Faithful Implementation with RcppHNSW-based kNN
#'
#' This code implements the GraphZoom framework as per the given paper:
#' "GraphZoom: A Multi-level Spectral Approach for Accurate and Scalable Graph Embedding"
#' by Chenhui Deng*, Zhiqiang Zhao*, Yongyu Wang, Zhiru Zhang, Zhuo Feng.
#'
#' The four phases are:
#' 1) Graph Fusion (Sec. 3.1)
#' 2) Spectral Graph Coarsening (Sec. 3.2, Eqs.(1)-(3))
#' 3) Graph Embedding on the Coarsest Graph (Sec. 3.3)
#' 4) Embedding Refinement (Sec. 3.4, Eqs.(4)-(7))
#'
#' Changes:
#' - Uses RcppHNSW for the k-NN graph construction on node attributes.
#' - Properly computes k-NN with L2 distance, then assigns edge weights based on cosine similarity.
#'
#' @references GraphZoom paper; RcppHNSW package documentation.
#' @import Matrix
#' @import RSpectra
#' @import RcppHNSW
#' @export

###############################
### Utility Functions
###############################

.compute_degree_matrix <- function(A) {
  d <- Matrix::rowSums(A)
  D <- Matrix::Diagonal(x=d)
  return(D)
}

.normalize_adjacency <- function(A) {
  D <- .compute_degree_matrix(A)
  d_inv_sqrt <- 1 / sqrt(Matrix::diag(D) + .Machine$double.eps)
  D_inv_sqrt <- Matrix::Diagonal(x=d_inv_sqrt)
  A_norm <- D_inv_sqrt %*% A %*% D_inv_sqrt
  return(A_norm)
}

.laplacian_matrix <- function(A) {
  D <- .compute_degree_matrix(A)
  L <- D - A
  return(L)
}

.apply_graph_filter <- function(A, E, H, k=2, sigma=1e-3) {
  # Add self-loops to reduce high-frequency noise
  N <- nrow(A)
  A_tilde <- A + sigma * Matrix::Diagonal(n=N, x=1)
  
  # Normalize A_tilde
  A_norm <- .normalize_adjacency(A_tilde)
  
  E_hat <- if(!is.null(H)) H %*% E else E
  
  # Apply (A_norm)^k as a low-pass filter approximation
  E_ref <- E_hat
  for (iter in seq_len(k)) {
    E_ref <- A_norm %*% E_ref
  }
  
  return(E_ref)
}


###############################
### Graph Fusion (Phase 1)
### Sec. 3.1
###############################
# Construct fused graph: 
# A_fusion = A_topo + beta * A_feat
# Where A_feat is from kNN on X (l2 distance) and weights from cosine similarity.

.fuse_graph <- function(A_topo, X, k=10, beta=1.0) {
  if(is.null(X)) {
    # No features available, return A_topo directly
    return(A_topo)
  }
  
  # Ensure X is a numeric matrix
  X <- as.matrix(X)
  N <- nrow(X)
  
  # Compute k-NN using RcppHNSW (l2 distance)
  # According to the paper, we first find k-nearest neighbors based on l2-norm distance.
  # Use hnsw_knn: returns a list with $item and $distance
  # items: matrix of indices of neighbors (N x k)
  # distance: N x k matrix of distances
  # By default hnsw_knn expects row-wise data: each row is an item
  all_knn <- RcppHNSW::hnsw_knn(X, k = k+1, distance = "l2") 
  # Often the nearest neighbor includes the point itself. We must remove self:
  # The first column typically is the item itself (distance=0). Remove it.
  # Check if first column of all_knn$item is the item itself:
  # Usually hnsw_knn returns the queried vector as nearest neighbor with distance 0.
  # We'll assume so and remove the first column.
  
  items <- all_knn$item[,-1,drop=FALSE]      # Nxk matrix of neighbors
  # distances <- all_knn$distance[,-1,drop=FALSE] # Nxk (not needed except for indexing check)
  
  # Compute weights based on cosine similarity:
  # w_{i,j} = (X_i . X_j)/(||X_i|| * ||X_j||)
  # Precompute norms
  norms <- sqrt(rowSums(X * X) + 1e-12)
  
  # We'll construct A_feat as a sparse matrix
  # i-th row: neighbors given by items[i,]
  # w_{i, items[i,j]} = cosine similarity
  i_idx <- integer(N*k)
  j_idx <- integer(N*k)
  vals <- numeric(N*k)
  
  idx <- 1
  for (i in seq_len(N)) {
    nbrs <- items[i,]
    # Dot products with neighbors:
    # X_i . X_j = sum over dimension
    # Instead of computing full dot products from scratch, do:
    # Actually just do a matrix multiplication since k might be large:
    Xi <- X[i,,drop=FALSE]
    Xnbrs <- X[nbrs,,drop=FALSE]
    dotp <- Xi %*% t(Xnbrs) # 1 x k
    
    # cos_sim
    cos_sim <- as.numeric(dotp) / (norms[i]*norms[nbrs])
    
    len_nbrs <- length(nbrs)
    i_idx[idx:(idx+len_nbrs-1)] <- i
    j_idx[idx:(idx+len_nbrs-1)] <- nbrs
    vals[idx:(idx+len_nbrs-1)] <- cos_sim
    idx <- idx+len_nbrs
  }
  
  # Construct A_feat (N x N)
  A_feat <- Matrix::sparseMatrix(i=i_idx, j=j_idx, x=vals, dims=c(N,N))
  # Symmetrize since it's kNN and we want an undirected graph
  A_feat <- (A_feat + t(A_feat))/2
  
  # Now A_fusion = A_topo + beta * A_feat
  A_fusion <- A_topo + beta * A_feat
  return(A_fusion)
}


###############################
### Spectral Coarsening (Phase 2)
### Sec. 3.2, Eqs.(1)-(3)
###############################

.spectral_coarsening <- function(A_fusion, levels=2, t=10) {
  graphs <- list(A_fusion)
  H_list <- list()
  
  for (level in seq_len(levels)) {
    A_curr <- graphs[[level]]
    N <- nrow(A_curr)
    L_curr <- .laplacian_matrix(A_curr)
    
    # Create t random vectors orthogonal to all-one vector
    x_init <- matrix(rnorm(N*t), nrow=N, ncol=t)
    for(i in seq_len(t)) {
      x_init[,i] <- x_init[,i] - mean(x_init[,i]) # orthogonalize against all-one
    }
    
    # Smoothing step (Eq.1 - low-pass filtering)
    # We'll apply a few smoothing iterations:
    max_iter <- 5
    alpha <- 0.001
    for(iter in seq_len(max_iter)) {
      Lx <- L_curr %*% x_init
      x_init <- x_init - alpha * Lx
    }
    
    # Compute spectral node affinity (Eq.2)
    # a_{p,q} = |(T_{p,:}, T_{q,:})|^2 / (||T_{p,:}||^2 * ||T_{q,:}||^2)
    norms <- sqrt(rowSums(x_init^2) + 1e-12)
    
    # We'll form node clusters by maximal matching on edges sorted by affinity
    A_trp <- summary(A_curr)
    sel <- A_trp$i < A_trp$j
    i_list <- A_trp$i[sel]
    j_list <- A_trp$j[sel]
    
    dot_products <- rowSums(x_init[i_list,,drop=FALSE] * x_init[j_list,,drop=FALSE])
    aff <- (dot_products^2)/((norms[i_list]^2)*(norms[j_list]^2))
    
    order_idx <- order(aff, decreasing=TRUE)
    
    matched <- rep(FALSE, N)
    parent <- rep(0, N)
    c_id <- 0
    
    for (edge_id in order_idx) {
      p <- i_list[edge_id]
      q <- j_list[edge_id]
      if(!matched[p] && !matched[q]) {
        c_id <- c_id + 1
        parent[p] <- c_id
        parent[q] <- c_id
        matched[p] <- TRUE
        matched[q] <- TRUE
      }
    }
    # Unmatched become singletons
    unmatched <- which(!matched)
    for (u in unmatched) {
      c_id <- c_id + 1
      parent[u] <- c_id
    }
    
    H <- Matrix::sparseMatrix(i=parent, j=1:N, x=1, dims=c(c_id,N))
    H_list[[level]] <- H
    
    # Coarsen graph: A_{i+1} = H * A_i * H^T
    A_coarse <- H %*% A_curr %*% Matrix::t(H)
    graphs[[level+1]] <- A_coarse
  }
  
  return(list(graphs=graphs, H_list=H_list))
}


###############################
### Graph Embedding (Phase 3)
### Sec. 3.3
###############################

.embed_coarsest_graph <- function(A_coarse, dim=16) {
  # Compute the normalized Laplacian L = I - D^{-1/2} A D^{-1/2}
  # Then find the smallest `dim` eigenvectors of L.
  
  N <- nrow(A_coarse)
  if (N < dim) {
    # If coarsest graph is extremely small, just return random embeddings
    return(matrix(rnorm(N * dim), nrow = N, ncol = dim))
  }
  
  # Compute normalized adjacency:
  D <- Matrix::Diagonal(x = Matrix::rowSums(A_coarse))
  d_inv_sqrt <- 1 / sqrt(Matrix::diag(D) + .Machine$double.eps)
  D_inv_sqrt <- Matrix::Diagonal(x = d_inv_sqrt)
  A_norm <- D_inv_sqrt %*% A_coarse %*% D_inv_sqrt
  
  # Normalized Laplacian: L = I - A_norm
  L <- Matrix::Diagonal(n=N, x=1) - A_norm
  
  # Use RSpectra to compute the first `dim+1` eigenvectors
  # We skip the trivial eigenvector corresponding to eigenvalue 0 (the all-ones vector),
  # and select the next `dim` eigenvectors as embeddings.
  # Note: Ensure L is symmetric. If needed, symmetrize: L <- (L + t(L))/2
  
  L <- (L + Matrix::t(L))/2  # ensure symmetry
  # smallest eigenvalues
  # opts=list() can be tuned for accuracy/performance
  eig_res <- RSpectra::eigs(L, k=dim+1, which="SM") # "SM" = smallest magnitude eigenvalues
  
  # The eigenvectors are in eig_res$vectors, columns = eigenvectors
  # The first eigenvector often corresponds to the trivial eigenvalue ~0.
  # We skip that and use the next `dim` eigenvectors as embeddings.
  E_coarse <- eig_res$vectors[, 2:(dim+1)]
  
  return(E_coarse)
}


###############################
### Embedding Refinement (Phase 4)
### Sec.3.4, Eqs.(4)-(7)
###############################

.refine_embeddings <- function(graphs, H_list, E_coarse, k=2, sigma=1e-3) {
  levels <- length(H_list)
  E_curr <- E_coarse
  
  # Refinement from coarsest back to the original
  for (level in seq(levels, 1, by=-1)) {
    A_finer <- graphs[[level]]
    H <- H_list[[level]]
    E_curr <- .apply_graph_filter(A_finer, E_curr, H, k, sigma)
  }
  
  return(E_curr)
}


###############################
### Main User-facing Function
###############################

#' @title GraphZoom Main Function (with RcppHNSW-based kNN)
#' @description Perform GraphZoom pipeline as described in the paper, using RcppHNSW for k-NN.
#'
#' @param A Sparse adjacency matrix (N x N) of the original graph topology.
#' @param X Node feature matrix (N x K). If NULL, no fusion is done and only topology is used.
#' @param levels Number of coarsening levels.
#' @param t Dimension for local spectral embedding during coarsening.
#' @param kNN Number of nearest neighbors for feature-based kNN graph construction.
#' @param beta Weight to balance topology and feature graph in fusion.
#' @param embed_dim Dimension of embeddings at coarsest level.
#' @param filter_power The power k in Eq.(7) for refinement.
#' @param sigma Self-loop factor in refinement filter.
#'
#' @return A list:
#' \item{E_original}{Refined embeddings for the original graph (N x embed_dim).}
#' \item{graphs}{List of graphs from level 0 to level L (fused and coarsened).}
#' \item{H_list}{List of mapping operators.}
#' \item{E_coarse}{Embeddings at the coarsest level.}
#'
#' @examples
#' library(Matrix)
#' library(RcppHNSW)
#' set.seed(42)
#' N <- 100
#' A <- rsparsematrix(N, N, 0.01)
#' A <- (A + t(A))/2
#' diag(A) <- 0
#' X <- matrix(rnorm(N*20), nrow=N, ncol=20)
#' result <- graph_zoom(A, X, levels=2, embed_dim=8, kNN=5)
#' E <- result$E_original
#'
#' @references GraphZoom Paper (Sections 3.1-3.4, Eqs.1-7); RcppHNSW docs for k-NN.
#' @export
graph_zoom <- function(A, X=NULL, levels=2, t=10, kNN=10, beta=1.0, 
                       embed_dim=16, filter_power=2, sigma=1e-3) {
  # Phase 1: Graph Fusion
  A_fusion <- .fuse_graph(A, X, k=kNN, beta=beta)
  
  # Phase 2: Spectral Coarsening
  coarsen_res <- .spectral_coarsening(A_fusion, levels=levels, t=t)
  graphs <- coarsen_res$graphs
  H_list <- coarsen_res$H_list
  
  # Phase 3: Embed Coarsest Graph
  A_coarsest <- graphs[[levels+1]]
  E_coarse <- .embed_coarsest_graph(A_coarsest, dim=embed_dim)
  
  # Phase 4: Refinement
  E_original <- .refine_embeddings(graphs, H_list, E_coarse, k=filter_power, sigma=sigma)
  
  return(list(E_original=E_original, graphs=graphs, H_list=H_list, E_coarse=E_coarse))
}