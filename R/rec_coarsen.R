#' Randomized Edge Contraction (REC) Graph Coarsening
#'
#' @description
#' Implements the randomized edge contraction (REC) coarsening method as described in the paper:
#' 
#' * A. Loukas and P. Vandergheynst, *"Spectrally Approximating Large Graphs with Smaller Graphs,"* 
#'   Proceedings of the 35th International Conference on Machine Learning (ICML), 2018.
#' 
#' Given a weighted graph represented by a sparse adjacency matrix, this function attempts to 
#' coarsen the graph by contracting edges at random according to a given edge potential. The result 
#' is a reduced (coarsened) graph that preserves the spectrum of the original Laplacian up to certain 
#' guarantees detailed in the paper.
#' 
#' @section Method:
#' The coarsening construction follows the framework detailed in Section 2 and Algorithm 1 of the paper. 
#' The main steps are:
#' 
#' 1. **Set-up**: Given a graph \eqn{G=(\mathcal{V},\mathcal{E},W)} with \eqn{N=|\mathcal{V}|}, 
#'    we form its Laplacian \eqn{L} defined in Equation (1) of the paper. 
#'    
#' 2. **Edge Potentials**: Assign a potential \eqn{\phi_{ij}} to each edge \eqn{e_{ij} \in \mathcal{E}}. 
#'    (See Section 3 of the paper). A common choice is the heavy-edge potential \eqn{\phi_{ij}=w_{ij}}.
#'    
#' 3. **Randomized Edge Contraction (REC)**: We perform up to \eqn{T} iterations (see Algorithm 1 in the paper):
#'    - At each iteration, from the set of candidate edges \eqn{\mathcal{C}}, 
#'      select one edge \eqn{e_{ij}} with probability proportional to its potential \eqn{\phi_{ij}}.
#'    - Contract \eqn{e_{ij}} (merge vertices \eqn{v_i} and \eqn{v_j}) to form a coarse vertex.
#'    - Remove all edges in the "edge neighborhood" \eqn{\mathcal{N}_{ij}} (those incident on \eqn{v_i} or \eqn{v_j}) from \eqn{\mathcal{C}}.
#'    - If no edge is selected at some iteration, we proceed to the next iteration (or stop if we have reached iteration \eqn{T}).
#'    
#' 4. **Coarsening Matrix**: Construct the coarsening matrix \eqn{C} as described in Section 2.1 (Equation (2) and discussion around it). 
#'    Each contracted set of original vertices forms one coarse vertex. Within each coarse vertex, the coarsening weights are constant 
#'    and normalized such that the norm of each block is 1.
#'
#' 5. **Result**: The output includes:
#'    - The coarsening matrix \eqn{C} (\eqn{n \times N}, with \eqn{n < N}),
#'    - The coarsened adjacency matrix \eqn{W_c},
#'    - The coarsened Laplacian \eqn{L_c = C L C^\top}.
#'
#' @section References:
#' 
#' - Equation (1): Definition of Laplacian.
#' - Equation (2): Definition of coarsened Laplacian.
#' - Algorithm 1: Randomized Edge Contraction.
#' - Sections 2 and 3: Construction and analysis of the coarsening.
#' 
#' @param W A sparse \eqn{N \times N} adjacency matrix of the original graph. Must be symmetric and non-negative.
#' @param T The number of iterations for randomized edge contraction. A larger \eqn{T} increases the chance of forming a large matching.
#' @param phi A vector or matrix of the same dimension as W (or length = number of edges) that encodes the edge potentials \eqn{\phi_{ij}}. 
#'   If \code{NULL}, defaults to heavy-edge potentials \eqn{\phi_{ij} = w_{ij}}.
#' @param seed An integer seed for reproducibility of random sampling.
#' @param deterministic_first_edge Logical, if TRUE selects the first edge in the candidate set for the first iteration.
#'   This is primarily for testing purposes. Default is FALSE.
#' @param verbose Logical, if TRUE prints debug information during the coarsening process.
#' @param use_cpp Logical. If TRUE, uses the C++ implementation for faster computation. Default is FALSE.
#' @param allow_multiple_passes Logical. If TRUE, allows multiple passes of coarsening. Default is FALSE.
#' 
#' @return A list with elements:
#' \describe{
#'   \item{C}{The \eqn{n \times N} coarsening matrix, where \eqn{n} is the number of coarse vertices.}
#'   \item{W_c}{The \eqn{n \times n} adjacency matrix of the coarsened graph.}
#'   \item{L_c}{The \eqn{n \times n} coarsened Laplacian matrix \eqn{L_c = C L C^\top}.}
#'   \item{mapping}{An integer vector of length \eqn{N} indicating the coarse vertex each original vertex belongs to.}
#' }
#'
#' @examples
#' # Example on a small random graph
#' set.seed(42)
#' library(Matrix)
#' N <- 10
#' # Create a random adjacency matrix for a weighted undirected graph
#' A <- rsparsematrix(N, N, density=0.2, symmetric=TRUE, rand.x=function(n) runif(n,0,1))
#' diag(A) <- 0
#' # Coarsen the graph
#' result <- rec_coarsen(A, T=10)
#' str(result)
#' 
#' @export
#' @importFrom Matrix Diagonal isSymmetric
rec_coarsen <- function(W, T=100, phi=NULL, seed=NULL, deterministic_first_edge=FALSE, 
                        verbose=FALSE, use_cpp=FALSE, allow_multiple_passes=FALSE) {
  if(!inherits(W, "dgCMatrix") && !inherits(W, "dgTMatrix") && !inherits(W, "dsCMatrix")) {
    W <- as(W, "dgCMatrix")
  }
  
  if(!is.null(seed)) set.seed(seed)
  
  N <- nrow(W)
  if(ncol(W)!=N) stop("W must be square.")
  
  # Ensure symmetry
  if(!isSymmetric(W)) stop("W must be symmetric.")
  
  if(!allow_multiple_passes) {
    if(use_cpp) {
      if(verbose) cat("Using C++ implementation...\n")
      result <- rec_coarsen_impl(W, T, phi, deterministic_first_edge, verbose, seed)
      if(verbose) cat("Algorithm stopped because:", result$stop_reason, "\n")
      result$stop_reason <- NULL
      return(result)
    }
    
    # Extract edges (i < j to avoid duplicates)
    # We'll store edges in a data.frame: (i,j,w,phi)
    # i,j are 1-based indices of vertices
    # Convert sparse matrix to a triplet form
    Wi <- W@i + 1L
    Wp <- W@p
    Wj <- integer(length(Wi))
    for (col in seq_len(N)) {
      start_idx <- Wp[col] + 1L
      end_idx <- Wp[col+1L]
      if (start_idx <= end_idx) {
        Wj[start_idx:end_idx] <- col
      }
    }
    
    # Edges with i<j only
    sel <- Wi < Wj
    e_i <- Wi[sel]
    e_j <- Wj[sel]
    w_e <- W@x[sel]
    
    edges <- data.frame(
      i = e_i,
      j = e_j,
      w = w_e,
      stringsAsFactors = FALSE
    )
    
    # If phi is NULL, use heavy-edge potentials phi_ij = w_ij
    # Otherwise, phi must correspond to each edge in edges
    if(is.null(phi)) {
      phi <- edges$w
    } else {
      # Check length
      if(length(phi) != nrow(edges)) stop("phi length must match the number of edges.")
    }
    edges$phi <- phi
    
    # Initially, each vertex is its own group
    parent <- seq_len(N) # union-find structure
    size <- rep(1L, N)
    
    # Union-Find helper functions
    find_set <- function(x) {
      if (parent[x] != x) {
        parent[x] <- find_set(parent[x])  # Path compression
      }
      parent[x]
    }
    
    union_set <- function(a, b) {
      a_root <- find_set(a)
      b_root <- find_set(b)
      
      if (a_root != b_root) {
        if (size[a_root] < size[b_root]) {
          # Swap to ensure a_root is the larger set
          tmp <- a_root
          a_root <- b_root
          b_root <- tmp
        }
        # Merge b into a
        parent[b_root] <<- a_root  # Use <<- to modify parent in parent environment
        size[a_root] <<- size[a_root] + size[b_root]  # Use <<- to modify size in parent environment
        if(verbose) {
          cat("After union: parent =", parent, "\n")
          cat("After union: size =", size, "\n")
        }
      }
    }
    
    # Candidate set C initially all edges
    cand <- rep(TRUE, nrow(edges))
    
    # Precompute neighborhood sets
    vertex_edge_map <- vector("list", N)
    for(eid in seq_len(nrow(edges))) {
      vertex_edge_map[[edges$i[eid]]] <- c(vertex_edge_map[[edges$i[eid]]], eid)
      vertex_edge_map[[edges$j[eid]]] <- c(vertex_edge_map[[edges$j[eid]]], eid)
    }
    
    total_phi <- sum(edges$phi[cand])
    
    # Track iterations used
    iterations_used <- 0
    
    # Perform up to T iterations
    for(iter in seq_len(T)) {
      iterations_used <- iter  # Track current iteration
      
      cands_idx <- which(cand)
      if(length(cands_idx) == 0) break
      
      chosen_edge_id <- NA_integer_
      
      if(deterministic_first_edge && iter == 1) {
        chosen_edge_id <- cands_idx[1]
        if(verbose) cat("Deterministically chose edge", chosen_edge_id, 
                       "connecting vertices", edges$i[chosen_edge_id], "and", edges$j[chosen_edge_id], "\n")
      } else {
        p <- edges$phi[cands_idx] / total_phi
        r <- runif(1)
        cum_p <- cumsum(p)
        
        if(r <= cum_p[length(cum_p)]) {
          chosen_edge_id <- cands_idx[which(r <= cum_p)[1]]
          if(verbose) cat("Probabilistically chose edge", chosen_edge_id, 
                         "connecting vertices", edges$i[chosen_edge_id], "and", edges$j[chosen_edge_id], "\n")
        } else {
          if(verbose) cat("No edge selected in iteration", iter, "\n")
          next
        }
      }
      
      # Contract chosen edge
      ei <- edges$i[chosen_edge_id]
      ej <- edges$j[chosen_edge_id]
      
      # Find leaders of sets
      pi <- find_set(ei)
      pj <- find_set(ej)
      
      if(verbose) {
        cat("Before merge: parent =", parent, "\n")
        cat("Before merge: size =", size, "\n")
      }
      
      if(pi != pj) {
        if(verbose) cat("Merging sets with leaders", pi, "and", pj, "\n")
        union_set(pi, pj)
        
        # Remove neighborhood edges from candidate set
        Ni <- vertex_edge_map[[ei]]
        Nj <- vertex_edge_map[[ej]]
        Nij <- unique(c(Ni, Nj))
        cand[Nij] <- FALSE
        
        if(verbose) cat("Removed edges", Nij, "from candidate set\n")
        
        total_phi <- sum(edges$phi[cand])
        if(total_phi == 0) break
      } else {
        if(verbose) cat("Vertices", ei, "and", ej, "already in same set\n")
      }
    }
    
    # After contraction, ensure all paths are compressed
    if(verbose) {
      cat("\nBefore final compression: parent =", parent, "\n")
      cat("Before final compression: size =", size, "\n")
    }
    
    # Compress paths and create mapping
    for(v in seq_len(N)) {
      parent[v] <- find_set(v)
    }
    
    if(verbose) {
      cat("After final compression: parent =", parent, "\n")
      cat("Final size array:", size, "\n")
    }
    
    # Get unique sets and create mapping
    unique_sets <- unique(parent)
    new_id_map <- seq_along(unique_sets)
    names(new_id_map) <- as.character(unique_sets)
    
    mapping <- integer(N)
    for(v in seq_len(N)) {
      mapping[v] <- new_id_map[as.character(parent[v])]
    }
    
    if(verbose) {
      cat("\nFinal mapping:", mapping, "\n")
      cat("Unique sets:", unique_sets, "\n")
      cat("Number of coarse vertices:", length(unique_sets), "\n")
    }
    
    n <- length(unique_sets) # number of coarse vertices
    
    # Construct coarsening matrix C
    # For each coarse vertex, find the vertices merged into it
    groups <- split(seq_len(N), mapping)
    
    # Each block c_j is 1/sqrt(n_j) for vertices in that coarse vertex, 0 elsewhere.
    # We'll construct C as a sparse matrix
    # C is n x N
    # The i-th row of C corresponds to coarse vertex i
    row_idx <- integer(N)
    col_idx <- integer(N)
    x_vals <- numeric(N)
    
    pos <- 1L
    for(ci in seq_along(groups)) {
      grp <- groups[[ci]]
      sz <- length(grp)
      val <- 1/sqrt(sz)
      row_idx[pos:(pos+sz-1)] <- ci
      col_idx[pos:(pos+sz-1)] <- grp
      x_vals[pos:(pos+sz-1)] <- val
      pos <- pos + sz
    }
    
    C <- sparseMatrix(i=row_idx[1:(pos-1)], j=col_idx[1:(pos-1)], 
                    x=x_vals[1:(pos-1)], dims=c(n,N))
    
    # Construct coarsened adjacency W_c = C W C^T
    CW <- C %*% W
    W_c <- CW %*% t(C)
    
    # Compute L_c = C L C^T
    d <- rowSums(W)
    L <- Diagonal(x=d) - W
    L_c <- C %*% L %*% t(C)
    
    list(
      C = C,
      W_c = W_c,
      L_c = L_c,
      mapping = mapping,
      iterations_used = iterations_used
    )
  } else {
    # Multiple passes logic
    remaining_T <- T
    current_W <- W
    current_result <- NULL
    pass <- 1
    
    while(remaining_T > 0) {
      if(verbose) cat(sprintf("\nStarting pass %d with T=%d\n", pass, remaining_T))
      
      # Run one pass
      pass_result <- if(use_cpp) {
        result <- rec_coarsen_impl(current_W, remaining_T, phi, 
                                 deterministic_first_edge, verbose, seed)
        if(verbose) cat("Pass", pass, "stopped because:", result$stop_reason, "\n")
        iterations_used <- result$iterations_used
        result$stop_reason <- NULL
        result$iterations_used <- NULL
        result
      } else {
        # Call the R implementation with remaining_T
        rec_coarsen(current_W, T=remaining_T, phi=phi, seed=seed,
                   deterministic_first_edge=deterministic_first_edge,
                   verbose=verbose, use_cpp=FALSE, allow_multiple_passes=FALSE)
      }
      
      # If no further coarsening possible or only one vertex remains, stop
      if(nrow(pass_result$W_c) == nrow(current_W) || nrow(pass_result$W_c) == 1) {
        if(verbose) cat("No further coarsening possible\n")
        break
      }
      
      # Compose with previous results if they exist
      if(is.null(current_result)) {
        current_result <- pass_result
      } else {
        current_result <- compose_coarsening(current_result, pass_result)
      }
      
      # Update for next pass
      current_W <- pass_result$W_c
      remaining_T <- remaining_T - iterations_used
      pass <- pass + 1
      
      if(verbose) {
        cat(sprintf("After pass %d: vertices reduced from %d to %d, remaining T: %d\n",
                  pass-1, nrow(W), nrow(current_W), remaining_T))
      }
    }
    
    return(current_result)
  }
}


#' Compose Two Coarsening Passes
#'
#' Given the results of two consecutive passes of the REC coarsening, this function composes
#' them to produce a single equivalent coarsening mapping and corresponding coarsened matrices.
#'
#' @param pass1 A list containing `C`, `W_c`, `L_c`, and `mapping` from the first coarsening pass.
#' @param pass2 A list containing `C`, `W_c`, `L_c`, and `mapping` from the second coarsening pass.
#'
#' @return A list with elements:
#' \describe{
#'   \item{C_final}{The final coarsening matrix mapping from the original graph to the final coarsened graph.}
#'   \item{W_c_final}{The final coarsened adjacency matrix.}
#'   \item{L_c_final}{The final coarsened Laplacian matrix.}
#'   \item{mapping_final}{An integer vector mapping original vertices directly to the final coarsened vertices.}
#' }
#'
#' @examples
#' # Suppose you have run:
#' # result1 <- rec_coarsen(W, T=..., ...)
#' # result2 <- rec_coarsen(result1$W_c, T=..., ...)
#' # Now compose:
#' # composed_result <- compose_coarsening(result1, result2)
#'
compose_coarsening <- function(pass1, pass2) {
  # Extract components from pass1
  C1 <- pass1$C
  W_c1 <- pass1$W_c
  L_c1 <- pass1$L_c
  mapping1 <- pass1$mapping
  
  # Extract components from pass2
  C2 <- pass2$C
  W_c2 <- pass2$W_c
  L_c2 <- pass2$L_c
  mapping2 <- pass2$mapping
  
  # Compose the coarsening matrices
  # Dimensions:
  # C1 is (n1 x N), C2 is (n2 x n1)
  # So C_final is (n2 x N)
  C_final <- C2 %*% C1
  
  # Compute W_c_final and L_c_final
  # Instead of recomputing from original W, we can use the relationships:
  # W_c_final = C2 * W_c1 * C2^T
  # L_c_final = C2 * L_c1 * C2^T
  W_c_final <- C2 %*% W_c1 %*% t(C2)
  L_c_final <- C2 %*% L_c1 %*% t(C2)
  
  # Compose the mappings:
  # mapping1: original vertices -> intermediate vertices
  # mapping2: intermediate vertices -> final vertices
  # mapping_final[i] = mapping2[mapping1[i]]
  
  mapping_final <- mapping2[mapping1]
  
  list(
    C = C_final,
    W_c = W_c_final,
    L_c = L_c_final,
    mapping = mapping_final
  )
}