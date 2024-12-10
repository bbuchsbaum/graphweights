#' GraphZoom: Faithful Implementation with RcppHNSW-based kNN
#'
#' This code implements the GraphZoom framework as per the given paper:
#' "GraphZoom: A Multi-level Spectral Approach for Accurate and Scalable Graph Embedding"
#'
#' Added Feature:
#' - You can now provide a precomputed feature graph (A_feat). If provided, it will be used directly.
#'   Otherwise, if X is provided, it will construct the feature graph from X. If neither is provided,
#'   it will rely solely on the topology graph.
#'
#' The four phases are:
#' 1) Graph Fusion (Sec. 3.1)
#' 2) Spectral Graph Coarsening (Sec. 3.2)
#' 3) Graph Embedding on the Coarsest Graph (Sec. 3.3)
#' 4) Embedding Refinement (Sec. 3.4)
#'
#' @references GraphZoom paper; RcppHNSW documentation
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
  N <- nrow(A)
  A_tilde <- A + sigma * Matrix::Diagonal(n=N, x=1)
  
  A_norm <- .normalize_adjacency(A_tilde)
  E_hat <- if(!is.null(H)) t(H) %*% E else E
  
  for (iter in seq_len(k)) {
    E_hat <- A_norm %*% E_hat
  }
  
  return(E_hat)
}


###############################
### Graph Fusion (Phase 1)
### Sec. 3.1
###############################
.fuse_graph <- function(A_topo, X, A_feat=NULL, k=10, beta=1.0, verbose=FALSE) {
  N <- nrow(A_topo)
  
  # If A_feat is provided directly, use it.
  if(!is.null(A_feat)) {
    if(verbose) cat("Feature graph provided directly. Using A_feat.\n")
    # Check dimensions
    if(nrow(A_feat) != N || ncol(A_feat) != N) {
      stop("Provided A_feat must be of dimension N x N, matching A_topo.")
    }
    A_fusion <- A_topo + beta * A_feat
    return(A_fusion)
  }
  
  # If no A_feat is provided, but X is given, construct A_feat from X
  if(!is.null(X)) {
    if(verbose) {
      cat("Computing k-NN graph with k =", k, "\n")
      cat("Input feature matrix dimensions:", dim(X), "\n")
    }
    
    X <- as.matrix(X)
    
    all_knn <- RcppHNSW::hnsw_knn(X, k = k+1, distance = "l2")
    if(verbose) {
      cat("k-NN computation complete\n")
      cat("Processing k-NN results...\n")
    }
    
    items <- all_knn$idx[,-1,drop=FALSE]
    norms <- sqrt(rowSums(X * X) + 1e-12)
    
    if(verbose) cat("Computing cosine similarities...\n")
    
    i_idx <- integer(N*k)
    j_idx <- integer(N*k)
    vals <- numeric(N*k)
    
    idx_ptr <- 1
    for (i in seq_len(N)) {
      nbrs <- items[i,]
      if(length(nbrs) == 0) {
        if(verbose) cat("Warning: No neighbors found for node", i, "\n")
        next
      }
      
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
    
    if(verbose) cat("Constructing feature graph...\n")
    A_feat <- Matrix::sparseMatrix(i=i_idx, j=j_idx, x=vals, dims=c(N,N))
    A_feat <- (A_feat + t(A_feat))/2
    
    if(verbose) cat("Computing fused graph...\n")
    A_fusion <- A_topo + beta * A_feat
    
    if(verbose) {
      cat("Graph fusion complete\n")
      cat("Final graph dimensions:", dim(A_fusion), "\n")
      cat("Number of non-zero entries:", length(A_fusion@x), "\n")
    }
    return(A_fusion)
  }
  
  # If neither A_feat nor X is provided, just return A_topo
  if(verbose) cat("No features or feature graph provided. Using topology only.\n")
  return(A_topo)
}


###############################
### Spectral Coarsening (Phase 2)
### Sec. 3.2
###############################
.spectral_coarsening <- function(A_fusion, levels=2, t=10, verbose=FALSE) {
  graphs <- list(A_fusion)
  H_list <- list()
  
  for (level in seq_len(levels)) {
    A_curr <- graphs[[level]]
    N <- nrow(A_curr)
    L_curr <- .laplacian_matrix(A_curr)
    
    x_init <- matrix(rnorm(N*t), nrow=N, ncol=t)
    for(i in seq_len(t)) {
      x_init[,i] <- x_init[,i] - mean(x_init[,i])
    }
    
    max_iter <- 5
    alpha <- 0.001
    for(iter in seq_len(max_iter)) {
      Lx <- L_curr %*% x_init
      x_init <- x_init - alpha * Lx
    }
    
    norms <- sqrt(rowSums(x_init^2) + 1e-12)
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
    
    unmatched <- which(!matched)
    for (u in unmatched) {
      c_id <- c_id + 1
      parent[u] <- c_id
    }
    
    H <- Matrix::sparseMatrix(i=parent, j=1:N, x=1, dims=c(c_id,N))
    H_list[[level]] <- H
    
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
  N <- nrow(A_coarse)
  if (N < dim) {
    return(matrix(rnorm(N * dim), nrow = N, ncol = dim))
  }
  
  D <- Matrix::Diagonal(x = Matrix::rowSums(A_coarse))
  d_inv_sqrt <- 1 / sqrt(Matrix::diag(D) + .Machine$double.eps)
  D_inv_sqrt <- Matrix::Diagonal(x = d_inv_sqrt)
  A_norm <- D_inv_sqrt %*% A_coarse %*% D_inv_sqrt
  
  L <- Matrix::Diagonal(n=N, x=1) - A_norm
  L <- (L + Matrix::t(L))/2
  eig_res <- RSpectra::eigs(L, k=dim+1, which="SM")
  
  E_coarse <- eig_res$vectors[, 2:(dim+1)]
  return(E_coarse)
}


###############################
### Embedding Refinement (Phase 4)
### Sec.3.4
###############################
.refine_embeddings <- function(graphs, H_list, E_coarse, k=2, sigma=1e-3) {
  levels <- length(H_list)
  E_curr <- E_coarse
  
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
#' @title GraphZoom Main Function with Optional Direct Feature Graph
#' @description Perform GraphZoom pipeline as described in the paper.
#'
#' If `A_feat` is provided, it will be used directly for graph fusion instead of computing
#' a feature graph from X. If both `A_feat` and `X` are provided, `A_feat` takes precedence.
#'
#' @param A Sparse adjacency matrix (N x N) of the original graph topology.
#' @param X Node feature matrix (N x K). If NULL, no fusion from features is done unless A_feat is provided.
#' @param A_feat Optional precomputed feature graph (N x N). If provided, used directly for fusion.
#' @param levels Number of coarsening levels.
#' @param t Dimension for local spectral embedding during coarsening.
#' @param kNN Number of nearest neighbors for feature-based kNN (ignored if A_feat is given).
#' @param beta Weight to balance topology and feature graph in fusion.
#' @param embed_dim Dimension of embeddings at coarsest level.
#' @param filter_power The power k in Eq.(7) for refinement.
#' @param sigma Self-loop factor in refinement filter.
#' @param verbose If TRUE, print progress messages.
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
#' # Run with computed feature graph:
#' result <- graph_zoom(A, X, levels=2, embed_dim=8, kNN=5)
#'
#' # If you have your own A_feat:
#' # Suppose A_feat is a NxN matrix you constructed from features externally:
#' # result <- graph_zoom(A, A_feat=my_feature_graph, levels=2, embed_dim=8)
#'
#' E <- result$E_original
#'
#' @references GraphZoom Paper (Sections 3.1-3.4)
#' @export
graph_zoom <- function(A, X=NULL, A_feat=NULL, levels=2, t=10, kNN=10, beta=1.0, 
                       embed_dim=16, filter_power=2, sigma=1e-3, verbose=FALSE) {
  if(verbose) cat("Starting GraphZoom pipeline...\n")
  
  # Phase 1: Graph Fusion
  if(verbose) {
    cat("Phase 1: Graph Fusion\n")
    cat("Input dimensions - A:", dim(A), "\n")
    cat("X:", if(!is.null(X)) dim(X) else "NULL", "\n")
    if(!is.null(A_feat)) cat("A_feat provided directly\n")
    cat("Parameters - kNN:", kNN, "beta:", beta, "\n")
  }
  A_fusion <- .fuse_graph(A, X, A_feat=A_feat, k=kNN, beta=beta, verbose=verbose)
  
  # Phase 2: Spectral Coarsening
  if(verbose) {
    cat("Phase 2: Spectral Coarsening\n")
    cat("Parameters - levels:", levels, "t:", t, "\n")
  }
  coarsen_res <- .spectral_coarsening(A_fusion, levels=levels, t=t, verbose=verbose)
  graphs <- coarsen_res$graphs
  H_list <- coarsen_res$H_list
  
  # Phase 3: Embedding the Coarsest Graph
  if(verbose) {
    cat("Phase 3: Embedding Coarsest Graph\n")
    cat("Coarsest graph dimensions:", dim(graphs[[levels+1]]), "\n")
    cat("Target embedding dimension:", embed_dim, "\n")
  }
  A_coarsest <- graphs[[levels+1]]
  E_coarse <- .embed_coarsest_graph(A_coarsest, dim=embed_dim)
  
  # Phase 4: Refinement
  if(verbose) {
    cat("Phase 4: Refinement\n")
    cat("Parameters - filter_power:", filter_power, "sigma:", sigma, "\n")
  }
  E_original <- .refine_embeddings(graphs, H_list, E_coarse, k=filter_power, sigma=sigma)
  
  if(verbose) cat("GraphZoom pipeline completed.\n")
  
  return(list(E_original=E_original, graphs=graphs, H_list=H_list, E_coarse=E_coarse))
}