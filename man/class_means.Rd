% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/all_generic.R
\name{class_means}
\alias{class_means}
\title{Class Means}
\usage{
class_means(x, i, ...)
}
\arguments{
\item{x}{An object.}

\item{i}{The index of the class for which the mean is to be computed.}

\item{...}{Additional arguments passed to specific methods.}
}
\value{
A matrix or data frame representing the means of each class, the structure of which depends on the input object's class.
}
\description{
A generic function to compute the mean of each class.
}
\examples{
# Create an example object
obj <- data.frame(x = c(1, 2, 3, 4, 5), y = c(2, 3, 5, 7, 11))

# Compute class means
means <- class_means(obj, 1)

}
