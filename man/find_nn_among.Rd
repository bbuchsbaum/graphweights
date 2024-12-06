% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/all_generic.R
\name{find_nn_among}
\alias{find_nn_among}
\title{Find nearest neighbors among a subset}
\usage{
find_nn_among(x, ...)
}
\arguments{
\item{x}{An object of class "nnsearcher" or "class_graph".}

\item{...}{Further arguments passed to or from other methods.}
}
\value{
A nearest neighbors result object.
}
\description{
Find nearest neighbors among a subset
}
\examples{
# Create an object of class "nnsearcher"
nn <- nnsearcher()

# Find nearest neighbors among a subset
result <- find_nn_among(nn)

}
