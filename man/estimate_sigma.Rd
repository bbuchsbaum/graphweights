% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/knn_weights.R
\name{estimate_sigma}
\alias{estimate_sigma}
\title{Estimate Optimal Bandwidth Parameter for Heat Kernel}
\usage{
estimate_sigma(X, prop = 0.25, nsamples = 500, normalized = FALSE)
}
\arguments{
\item{X}{Numeric matrix where rows are instances and columns are features}

\item{prop}{Numeric value in (0,1) specifying the quantile to use (default: 0.25)}

\item{nsamples}{Integer specifying maximum number of samples to use (default: 500)}

\item{normalized}{Logical indicating whether distances should be computed on
normalized data (default: FALSE)}
}
\value{
Numeric value representing the estimated optimal bandwidth parameter
}
\description{
Estimates an appropriate bandwidth parameter (sigma) for the heat kernel
based on the distribution of pairwise distances in the data. This is a crucial
parameter that affects the scale at which similarities are measured.
}
\details{
The function estimates sigma by:
1. Computing pairwise distances between points (or a sample of points)
2. Finding the specified quantile of non-zero distances
This approach helps choose a sigma value that's appropriate for the scale
of the data.
}
\examples{
# Generate example data
set.seed(123)
X <- matrix(rnorm(1000 * 10), nrow = 1000)

# Estimate sigma with default parameters
sigma1 <- estimate_sigma(X)

# Estimate sigma with different quantile
sigma2 <- estimate_sigma(X, prop = 0.1)

# Estimate sigma for normalized data
sigma3 <- estimate_sigma(X, normalized = TRUE)

}
