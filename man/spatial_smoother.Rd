% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/spatial_weights.R
\name{spatial_smoother}
\alias{spatial_smoother}
\title{spatial_smoother}
\usage{
spatial_smoother(
  coord_mat,
  sigma = 5,
  nnk = 3^(ncol(coord_mat)),
  stochastic = TRUE
)
}
\arguments{
\item{coord_mat}{the matrix of coordinates}

\item{sigma}{the standard deviation of the smoother}

\item{nnk}{the number of neighbors in the kernel}

\item{stochastic}{make doubly stochastic}
}
\description{
spatial_smoother
}
\examples{
coord_mat <- as.matrix(expand.grid(x=1:10, y=1:10))
sm <- spatial_smoother(coord_mat, sigma=3)
sm2 <- spatial_smoother(coord_mat, sigma=6, stochastic=FALSE,nnk=50)
}
