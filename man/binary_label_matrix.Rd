% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/label_sim.R
\name{binary_label_matrix}
\alias{binary_label_matrix}
\title{Create a Binary Label Adjacency Matrix}
\usage{
binary_label_matrix(a, b = NULL, type = c("s", "d"))
}
\arguments{
\item{a}{A vector of labels for the first set of data points.}

\item{b}{A vector of labels for the second set of data points (default: NULL). If NULL, `b` will be set to `a`.}

\item{type}{A character specifying the type of adjacency matrix to create, either "s" for same labels or "d" for different labels (default: "s").}
}
\value{
A binary adjacency matrix with edges determined by the label relationships between `a` and `b`.
}
\description{
Constructs a binary adjacency matrix based on two sets of labels `a` and `b`, creating edges when `a` and `b` have the same or different labels, depending on the `type` parameter.
}
\examples{
data(iris)
a <- iris[,5]
bl <- binary_label_matrix(a, type="d")

}
