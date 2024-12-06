% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/class_graph.R
\name{discriminating_similarity}
\alias{discriminating_similarity}
\title{Compute Similarity Graph Weighted by Class Structure}
\usage{
discriminating_similarity(X, k = length(labels)/2, sigma, cg, threshold = 0.01)
}
\arguments{
\item{X}{A numeric matrix or data frame where rows are observations and columns are features.}

\item{k}{Integer; number of nearest neighbors to consider. Defaults to half the number of labels.}

\item{sigma}{Numeric; scaling factor for the heat kernel.}

\item{cg}{A class_graph object containing class structure information.}

\item{threshold}{Numeric; threshold value for the class graph (default: 0.01).}
}
\value{
A numeric matrix representing the weighted similarity graph.
}
\description{
Computes a similarity graph that incorporates class structure information.
This method is particularly effective for pattern recognition tasks such as
face recognition and handwriting digit classification.
}
\examples{
# Create example data
set.seed(123)
X <- matrix(rnorm(100 * 4), 100, 4)
labels <- factor(rep(1:5, each = 20))
cg <- class_graph(labels)

# Compute discriminating similarity
W <- discriminating_similarity(X, k = 10, sigma = 0.7, cg = cg)

}
\references{
Zhang, W., Xue, X., Sun, Z., Guo, Y., & Lu, H. (2007).
Local similarity and diversity preserving discriminant projection for face
and handwriting digits recognition.
Neurocomputing, 70(7-9), 1454-1461.
}
