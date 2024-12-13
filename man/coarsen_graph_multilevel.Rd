% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/multi_coarsen.R
\name{coarsen_graph_multilevel}
\alias{coarsen_graph_multilevel}
\title{Multi-level Graph Coarsening by Local Variation (Edge-based)}
\usage{
coarsen_graph_multilevel(
  A,
  n,
  k = 2,
  epsilon_prime = 0.5,
  max_levels = NULL,
  verbose = FALSE
)
}
\arguments{
\item{A}{A sparse adjacency matrix (dgCMatrix) of the undirected graph.}

\item{n}{Target number of vertices after coarsening (must be < N).}

\item{k}{Dimension of the subspace to preserve (e.g., first k eigenvectors of L).}

\item{epsilon_prime}{Error threshold for restricted spectral approximation (used as a heuristic).}

\item{max_levels}{Maximum levels allowed (default: log2(N/n)), just to prevent infinite loops.}

\item{verbose}{If TRUE, prints progress messages.}
}
\value{
A list:
\item{A_coarse}{The coarsened adjacency matrix.}
\item{P}{The overall reduction matrix mapping original N to final n.}
\item{info}{A list with final_epsilon, final_n, levels.}
}
\description{
Implements the multi-level edge-based local variation coarsening 
based on Loukas (2019). The method updates the Laplacian and a subspace U at each level,
ensuring stable re-orthonormalization and forcing at least one contraction if needed.
}
\references{
Loukas, A. (2019). "Graph Reduction with Spectral and Cut Guarantees."
Journal of Machine Learning Research, 20(116):1–42.
}
