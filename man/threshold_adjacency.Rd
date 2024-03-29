% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/knn_weights.R
\name{threshold_adjacency}
\alias{threshold_adjacency}
\title{Threshold Adjacency}
\usage{
threshold_adjacency(A, k = 5, type = c("normal", "mutual"), ncores = 1)
}
\arguments{
\item{A}{An adjacency matrix representing the graph.}

\item{k}{An integer specifying the number of neighbors to consider (default: 5).}

\item{type}{A character string indicating the type of k-nearest neighbors graph to compute. One of "normal" or "mutual" (default: "normal").}

\item{ncores}{An integer specifying the number of cores to use for parallel computation (default: 1).}
}
\value{
A sparse adjacency matrix containing only the specified number of nearest neighbors.
}
\description{
This function extracts the k-nearest neighbors from an existing adjacency matrix.
It returns a new adjacency matrix containing only the specified number of nearest neighbors.
}
\examples{
A <- matrix(runif(100), 10, 10)
A_thresholded <- threshold_adjacency(A, k = 5)

}
