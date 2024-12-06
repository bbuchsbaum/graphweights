% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/spatial_weights.R
\name{cross_weighted_spatial_adjacency}
\alias{cross_weighted_spatial_adjacency}
\title{Compute the Adjacency between Two Adjacency Matrices Weighted by Two Corresponding Feature Matrices.}
\usage{
cross_weighted_spatial_adjacency(
  coord_mat1,
  coord_mat2,
  feature_mat1,
  feature_mat2,
  wsigma = 0.73,
  alpha = 0.5,
  nnk = 27,
  maxk = nnk,
  weight_mode = c("binary", "heat"),
  sigma = 1,
  dthresh = sigma * 2.5,
  normalized = TRUE
)
}
\arguments{
\item{coord_mat1}{the first coordinate matrix (the query)}

\item{coord_mat2}{the second coordinate matrix (the reference)}

\item{feature_mat1}{the first feature matrix}

\item{feature_mat2}{the second feature matrix}

\item{wsigma}{the sigma for the feature heat kernel}

\item{alpha}{the mixing weight for the spatial distance (1=all spatial weighting, 0=all feature weighting)}

\item{nnk}{the maximum number of spatial nearest neighbors to include}

\item{maxk}{the maximum number of neighbors to include within spatial window}

\item{dthresh}{the threshold for the spatial distance}

\item{normalized}{whether to normalize the rows to sum to 1}
}
\description{
Computes the adjacency between two adjacency matrices weighted by two corresponding feature matrices.
}
\details{
The function first computes the spatial adjacency matrix between the two sets of points,
and then weights the matrix using the feature matrices.
}
\examples{

coords <- as.matrix(expand.grid(1:5, 1:5))
fmat1 <- matrix(rnorm(5*25), 25, 5)
fmat2 <- matrix(rnorm(5*25), 25, 5)

adj <- cross_weighted_spatial_adjacency(coords, coords, fmat1, fmat2)
}
