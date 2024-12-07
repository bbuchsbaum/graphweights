% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/label_sim.R
\name{create_label_matrix}
\alias{create_label_matrix}
\title{Create a Label Similarity Matrix}
\usage{
create_label_matrix(
  a,
  b = NULL,
  type = c("sparse", "dense", "custom"),
  similarity_fn = NULL,
  threshold = NULL,
  mode = c("similarity", "distance")
)
}
\arguments{
\item{a}{Vector of labels for the first set of data points}

\item{b}{Optional vector of labels for the second set. If NULL, uses 'a'}

\item{type}{Type of relationship to compute:
* "sparse": Memory-efficient implementation (formerly binary_label_matrix)
* "dense": Fast implementation for smaller datasets (formerly label_matrix2)
* "custom": Flexible implementation with custom similarity (formerly label_matrix)}

\item{similarity_fn}{Optional custom similarity function}

\item{threshold}{Threshold for similarity values (default: NULL)}

\item{mode}{Character, either "similarity" or "distance" (default: "similarity")}
}
\value{
A sparse Matrix object containing the computed similarities/distances
}
\description{
Creates a similarity or dissimilarity matrix between categorical labels using the most
appropriate implementation based on the input characteristics and requirements.
}
\details{
This function provides a unified interface for creating label similarity matrices,
automatically selecting the most appropriate implementation based on:
* Input size and sparsity
* Memory constraints
* Custom similarity requirements

Available implementations:
* Sparse: Memory-efficient for large, sparse label relationships
* Dense: Fast for smaller datasets with many relationships
* Custom: Flexible implementation supporting arbitrary similarity functions
}
\examples{
# Basic sparse similarity matrix
labels <- factor(c("A", "B", "A", "C"))
sim_mat <- create_label_matrix(labels, type = "sparse")

# Dense matrix with custom similarity
custom_sim <- function(x, y) as.numeric(x == y) * 2
sim_mat2 <- create_label_matrix(labels, type = "dense",
                               similarity_fn = custom_sim)

}
