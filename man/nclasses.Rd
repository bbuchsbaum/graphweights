% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/all_generic.R
\name{nclasses}
\alias{nclasses}
\title{Number of Classes}
\usage{
nclasses(x)
}
\arguments{
\item{x}{An object.}
}
\value{
The number of classes in the input object.
}
\description{
A generic function to compute the number of classes in a graph.
}
\examples{
# Create an example graph
g <- neighbor_graph(make_ring(5))

# Compute number of classes
n_classes <- nclasses(g)

}
