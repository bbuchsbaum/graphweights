% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/neighbor_graph.R
\name{laplacian.neighbor_graph}
\alias{laplacian.neighbor_graph}
\title{Laplacian Matrix for neighbor_graph Objects}
\usage{
\method{laplacian}{neighbor_graph}(x, normalized = FALSE)
}
\arguments{
\item{x}{A neighbor_graph object}

\item{normalized}{Logical; whether to compute the normalized Laplacian
(default: FALSE)}
}
\value{
A sparse Matrix object representing the Laplacian matrix
}
\description{
Computes the Laplacian matrix of a neighbor_graph object, with optional
normalization.
}
\details{
The Laplacian matrix L is computed as L = D - A, where:
* D is the degree matrix (diagonal matrix of vertex degrees)
* A is the adjacency matrix

When normalized=TRUE, the normalized Laplacian L_norm = D^(-1/2) L D^(-1/2)
is computed instead.

The Laplacian matrix is fundamental in spectral graph theory and is used for:
* Spectral clustering
* Graph partitioning
* Dimensionality reduction
* Computing graph properties
}
\examples{
\dontrun{
# Create a neighbor graph
g <- make_ring(5)
ng <- neighbor_graph(g)

# Compute standard Laplacian
L <- laplacian(ng)

# Compute normalized Laplacian
L_norm <- laplacian(ng, normalized = TRUE)
}

}
\seealso{
* \code{\link{adjacency.neighbor_graph}} for the adjacency matrix
}
