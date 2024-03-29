% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/spatial_weights.R
\name{pairwise_adjacency}
\alias{pairwise_adjacency}
\title{Compute a pairwise adjacency matrix for multiple graphs}
\usage{
pairwise_adjacency(Xcoords, Xfeats, fself, fbetween)
}
\arguments{
\item{Xcoords}{A list of numeric matrices or data.frames containing the spatial coordinates of the nodes of each graph.}

\item{Xfeats}{A list of numeric matrices or data.frames containing the feature vectors for the nodes of each graph.}

\item{fself}{A function that computes similarity for nodes within the same graph (e.g., Xi_1, Xi_2).}

\item{fbetween}{A function that computes similarity for nodes across graphs (e.g., Xi_1, Xj_1).}
}
\value{
A sparse matrix representing the pairwise adjacency matrix for the input graphs.
}
\description{
This function computes a pairwise adjacency matrix for multiple graphs with a given set of spatial coordinates
and feature vectors. The function takes two user-defined functions to compute within-graph and between-graph
similarity measures.
}
