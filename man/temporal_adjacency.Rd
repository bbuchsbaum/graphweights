% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/temporal_weights.R
\name{temporal_adjacency}
\alias{temporal_adjacency}
\title{Compute Temporal Adjacency Matrix}
\usage{
temporal_adjacency(
  time,
  weight_mode = c("heat", "binary"),
  sigma = 1,
  window = 2
)
}
\arguments{
\item{time}{A numeric vector representing time points. Values should be ordered
and typically represent equally spaced time intervals.}

\item{weight_mode}{Character string specifying the weighting scheme:
* "heat": Gaussian heat kernel weights (exponential decay with distance)
* "binary": Binary weights (1 for points within window, 0 otherwise)
Default is "heat".}

\item{sigma}{Numeric value controlling the width of the Gaussian kernel when
weight_mode = "heat". Larger values create broader temporal connections.
Only used when weight_mode = "heat". Default is 1.}

\item{window}{Integer specifying the size of the temporal window. Determines
how many adjacent time points are connected. Must be >= 1. Default is 2.}
}
\value{
A sparse symmetric Matrix object representing temporal adjacency.
  The matrix has dimensions length(time) x length(time), where element (i,j)
  represents the connection strength between time points i and j.
}
\description{
Constructs a temporal adjacency matrix for a time series using a sliding window
approach. The adjacency can be binary or weighted using a heat kernel, capturing
the temporal relationships between time points.
}
\details{
The function creates an adjacency matrix where connections are established
between time points within a specified window. The strength of connections
can be either:
* Binary (1 for connected points, 0 otherwise)
* Heat kernel weighted (exponentially decaying with temporal distance)

The resulting matrix is symmetric and includes self-connections (diagonal
elements are 1). This structure is useful for temporal clustering, smoothing,
and pattern detection in time series data.
}
\examples{
\dontrun{
# Create evenly spaced time points
time <- seq(0, 10, by = 0.5)

# Compute binary temporal adjacency
adj_binary <- temporal_adjacency(time, 
                                weight_mode = "binary",
                                window = 2)

# Compute heat kernel weighted temporal adjacency
adj_heat <- temporal_adjacency(time,
                              weight_mode = "heat",
                              sigma = 1.5,
                              window = 3)
}

}
\seealso{
* \code{\link{temporal_autocor}} for computing temporal autocorrelation
* \code{\link{temporal_laplacian}} for computing temporal Laplacian matrices
}
