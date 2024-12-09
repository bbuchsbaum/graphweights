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
rec_coarsen <- function(W, T=100, phi=NULL, seed=NULL) {
  if(!inherits(W, "dgCMatrix") && !inherits(W, "dgTMatrix") && !inherits(W, "dsCMatrix")) {
    W <- as(W, "dgCMatrix")
  }
  
  if(!is.null(seed)) set.seed(seed)
  
  N <- nrow(W)
  if(ncol(W)!=N) stop("W must be square.")
  
  # Ensure symmetry
  if(!isSymmetric(W)) stop("W must be symmetric.")
  
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
  sel <- Wi < (Wj+1L)
  e_i <- Wi[sel]
  e_j <- (Wj+1L)[sel]
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
  # We'll perform edge contractions by merging groups
  parent <- seq_len(N) # union-find structure
  size <- rep(1L, N)
  
  # Union-Find helper functions
  find_set <- function(x) {
    while(parent[x] != x) {
      parent[x] <- parent[parent[x]]
      x <- parent[x]
    }
    x
  }
  
  union_set <- function(a, b) {
    if(size[a] < size[b]) {
      tmp <- a; a <- b; b <- tmp
    }
    parent[b] <- a
    size[a] <- size[a] + size[b]
  }
  
  # Candidate set C initially all edges
  # We'll store a logical vector indicating which edges are still candidate
  cand <- rep(TRUE, nrow(edges))
  
  # Precompute neighborhood sets: for each edge, which edges share a vertex?
  # This can be large, so we do it efficiently.
  # We'll create adjacency lists for edges indexed by vertex
  # vertex_edge_map[[v]] gives the indices of edges incident on vertex v
  vertex_edge_map <- vector("list", N)
  for(eid in seq_len(nrow(edges))) {
    vertex_edge_map[[edges$i[eid]]] <- c(vertex_edge_map[[edges$i[eid]]], eid)
    vertex_edge_map[[edges$j[eid]]] <- c(vertex_edge_map[[edges$j[eid]]], eid)
  }
  
  # total_phi = sum of phi over candidate edges
  # We'll update total_phi at each iteration if an edge is contracted
  total_phi <- sum(edges$phi[cand])
  
  # Perform up to T iterations
  for(iter in seq_len(T)) {
    cands_idx <- which(cand)
    if(length(cands_idx) == 0) break
    
    # Probability of picking each candidate edge
    p <- edges$phi[cands_idx] / total_phi
    # We may also continue with no edge selection
    # Probability that no edge is selected = 1 - sum(p)
    # Draw a random number
    r <- runif(1)
    cum_p <- cumsum(p)
    
    chosen_edge_id <- NA_integer_
    if(r <= cum_p[length(cum_p)]) {
      # An edge is selected
      chosen_edge_id <- cands_idx[which(r <= cum_p)[1]]
    } else {
      # continue without contraction
      # no edge selected this iteration
      next
    }
    
    # Contract chosen edge
    ei <- edges$i[chosen_edge_id]
    ej <- edges$j[chosen_edge_id]
    
    # Find leaders of sets
    pi <- find_set(ei)
    pj <- find_set(ej)
    if(pi != pj) {
      # Perform union (merge) of these sets
      union_set(pi, pj)
      
      # Remove neighborhood edges from candidate set
      # Neighborhood edges are all edges incident on ei or ej
      # Actually, we must consider the representative sets of ei and ej before merging
      # but since we just merged them, we consider old references
      Ni <- vertex_edge_map[[ei]]
      Nj <- vertex_edge_map[[ej]]
      Nij <- unique(c(Ni, Nj))
      
      # Mark these edges as not candidate
      for(ne in Nij) {
        if(cand[ne]) {
          cand[ne] <- FALSE
        }
      }
      
      # Update total_phi
      total_phi <- sum(edges$phi[cand])
      if(total_phi == 0) break
    }
    # else if pi == pj, the edge chosen connects vertices already in same set
    # effectively a self-loop in the coarse representation, we can ignore
  }
  
  # After contraction, we have a partition of the original vertices
  # Let's assign new IDs to each coarse vertex
  for(v in seq_len(N)) find_set(v) # path compression
  unique_sets <- unique(parent)
  new_id_map <- seq_along(unique_sets)
  names(new_id_map) <- as.character(unique_sets)
  
  mapping <- parent
  for(v in seq_len(N)) {
    mapping[v] <- new_id_map[as.character(parent[v])]
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
  
  C <- sparseMatrix(i=row_idx, j=col_idx, x=x_vals, dims=c(n,N))
  
  # Construct coarsened adjacency W_c
  # W_c = C W C^T
  # We could compute W_c directly by merging groups:
  # sum of edge weights between sets forms W_c.
  # An alternative: directly compute W_c from W and mapping.
  
  # Let's do a group-based approach:
  # The weight between two coarse vertices a,b is sum of w_ij for i in group a, j in group b.
  # We'll use a hash-based approach
  # NOTE: Since we often want L_c = C L C^T, 
  # W_c can also be computed as: W_c = C W C^T.
  # Direct multiplication might be expensive for large graphs, but let's do it straightforwardly.
  
  # Compute W_c = C W C^T
  # C W is n x N, then (C W) C^T is n x n
  CW <- C %*% W
  W_c <- CW %*% t(C)
  
  # Compute L_c = C L C^T
  # L = D - W, where D = diag(rowSums(W))
  d <- rowSums(W)
  L <- Diagonal(x=d) - W
  L_c <- C %*% L %*% t(C)
  
  list(
    C = C,
    W_c = W_c,
    L_c = L_c,
    mapping = mapping
  )
}