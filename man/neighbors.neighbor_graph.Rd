% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/neighbor_graph.R
\name{neighbors.neighbor_graph}
\alias{neighbors.neighbor_graph}
\title{Neighbors for neighbor_graph Objects}
\usage{
\method{neighbors}{neighbor_graph}(x, i)
}
\arguments{
\item{x}{A neighbor_graph object}

\item{i}{Optional; integer specifying the vertex index for which to find neighbors.
If missing, returns neighbors for all vertices}
}
\value{
If i is specified:
  * A vector of vertex indices representing neighbors of vertex i
If i is missing:
  * A list where each element contains the neighbors of the corresponding vertex
}
\description{
Retrieves the neighbors of a specific vertex or all vertices in a neighbor_graph
object.
}
\details{
For a given vertex i, returns all vertices that share an edge with i.
If no vertex is specified, returns a list of neighbors for all vertices
in the graph.
}
\examples{
\dontrun{
# Create a neighbor graph
g <- make_ring(5)
ng <- neighbor_graph(g)

# Get neighbors of vertex 1
n1 <- neighbors(ng, 1)

# Get all neighbors
all_n <- neighbors(ng)
}

}
\seealso{
* \code{\link{non_neighbors.neighbor_graph}} for finding non-neighboring vertices
* \code{\link{edges.neighbor_graph}} for getting all edges
}
