#' Construct a class_graph
#'
#' A graph in which members of the same class have edges.
#'
#' @param labels A vector of class labels.
#' @param sparse A logical value, indicating whether to use sparse matrices in the computation. Default is TRUE.
#'
#' @return A class_graph object, which is a list containing the following components:
#'   - adjacency: A matrix representing the adjacency of the graph.
#'   - params: A list of parameters used in the construction of the graph.
#'   - labels: A vector of class labels.
#'   - class_indices: A list of vectors, each containing the indices of elements belonging to a specific class.
#'   - class_freq: A table of frequencies for each class.
#'   - levels: A vector of unique class labels.
#'   - classes: A character string indicating the type of graph ("class_graph").
#'
#' @importFrom Matrix sparseVector tcrossprod Matrix t
#'
#' @examples
#' data(iris)
#' labels <- iris[,5]
#' cg <- class_graph(labels)
#'
#' @export
class_graph <- function(labels, sparse=TRUE) {
  labels <- as.factor(labels)
  out <- Reduce("+", lapply(levels(labels), function(lev) {
    kronecker(Matrix(labels==lev, sparse=sparse), t(Matrix(labels==lev, sparse=sparse)))
  }))

  ret <- neighbor_graph(
    out,
    params = list(weight_mode="binary", neighbor_mode="supervised"),
    labels=labels,
    class_indices=split(1:length(labels), labels),
    class_freq=table(labels),
    levels=unique(labels),
    classes="class_graph"
  )

  ret
}



#' @export
nclasses.class_graph <- function(x) {
  length(x$levels)
}


#' @export
class_means.class_graph <- function(x, X) {
  ret <- do.call(rbind, lapply(x$class_indices, function(i) {
    colMeans(X[i,,drop=FALSE])
  }))

  row.names(ret) <- names(x$class_indices)
  ret
}


heterogeneous_neighbors.class_graph <- function(x, X, k, weight_mode="binary", sigma=.7) {
  allind <- 1:nrow(X)
  cind <- x$class_indices
  ## this could be done set-wise, e.g. levels by level
  out <- do.call(rbind, lapply(cind, function(ci) {
    query <- X[ci,]
    others <- which(!(allind %in% ci))

    m <- weighted_knnx(X[others,,drop=FALSE], query,k=k,
                  FUN=get_neighbor_fun(weight_mode, len=nrow(X), sigma=sigma), as="index_sim")
    ## HER
  }))

  ng <- neighbor_graph(igraph::graph_from_data_frame(out,directed=FALSE))

}

homogeneous_neighbors.class_graph <- function(x, X, k, weight_mode="heat", sigma=1) {

  cind <- x$class_indices
  out <- do.call(rbind, lapply(cind, function(ci) {
    m <- weighted_knn(X[ci,,drop=FALSE], k=k,
                       FUN=get_neighbor_fun(weight_mode, len=length(cind), sigma=sigma), as="sparse")


    Tm <- as_triplet(m)
    Tm[,1] <- ci[Tm[,1]]
    Tm[,2] <- ci[Tm[,2]]
    Tm
  }))

  ng <- neighbor_graph(igraph::graph_from_data_frame(out,directed=FALSE))

}


#' Within-Class Neighbors for class_graph Objects
#'
#' Compute the within-class neighbors of a class_graph object.
#'
#' @param x A class_graph object.
#' @param ng A neighbor graph object.
#'
#' @return A neighbor_graph object representing the within-class neighbors of the input class_graph.
#'
#' @inheritParams graph_weights
#' @export
within_class_neighbors.class_graph <- function(x, ng) {
  Ac <- adjacency(x)
  An <- adjacency(ng)
  Aout <- Ac * An
  neighbor_graph(Aout)
}


#' @inheritParams graph_weights
#' @export
between_class_neighbors.class_graph <- function(x, ng) {
  Ac <- !adjacency(x)
  An <- adjacency(ng)
  Aout <- Ac * An
  neighbor_graph(Aout)
}


#' Compute Discriminating Distance for Similarity Graph
#'
#' This function computes a discriminating distance matrix for the similarity graph based on the class labels.
#' It adjusts the similarity graph by modifying the weights within and between classes, making it more suitable for
#' tasks like classification and clustering.
#'
#' @param X A numeric matrix or data frame containing the data points.
#' @param k An integer representing the number of nearest neighbors to consider. Default is the number of unique labels divided by 2.
#' @param sigma A numeric value representing the scaling factor for the heat kernel. If not provided, it will be estimated.
#' @param labels A factor or numeric vector containing the class labels for each data point.
#'
#' @return A discriminating distance matrix in the form of a numeric matrix.
#'
#' @examples
#' X <- matrix(rnorm(100*100), 100, 100)
#' labels <- factor(rep(1:5, each=20))
#' sigma <- 0.7
#' D <- discriminating_distance(X, k=length(labels)/2, sigma, labels)
#'
#' @export
discriminating_distance <- function(X, k=length(labels)/2, sigma,labels) {
  #Wknn <- graph_weights(X)

  if (missing(sigma)) {
    sigma <- estimate_sigma(X, prop=.1)
  }

  Wall <- graph_weights(X, k=k, weight_mode="euclidean", neighbor_mode="knn")

  Ww <- label_matrix2(labels, labels)
  Wb <- label_matrix2(labels, labels, type="d")

  Ww2 <- Wall * Ww
  Wb2 <- Wall * Wb

  wind <- which(Ww2 >0)
  bind <- which(Wb2 >0)

  hw <- inverse_heat_kernel(Wall[wind], sigma)
  hb <- inverse_heat_kernel(Wall[bind], sigma)

  Wall[wind] <- hw * (1-hw)
  Wall[bind] <- hb * (1+hb)
  Wall
}

#' Compute Similarity Graph Weighted by Class Structure
#'
#' This function computes a similarity graph that is weighted by the class structure of the data.
#' It is useful for preserving the local similarity and diversity within the data, making it
#' suitable for tasks like face and handwriting digits recognition.
#'
#' @param X A numeric matrix or data frame containing the data points.
#' @param k An integer representing the number of nearest neighbors to consider. Default is the number of unique labels divided by 2.
#' @param sigma A numeric value representing the scaling factor for the heat kernel.
#' @param cg A class_graph object computed from the labels.
#' @param threshold A numeric value representing the threshold for the class graph. Default is 0.01.
#'
#' @return A weighted similarity graph in the form of a matrix.
#'
#' @examples
#' X <- matrix(rnorm(100*100), 100, 100)
#' labels <- factor(rep(1:5, each=20))
#' cg <- class_graph(labels)
#' sigma <- 0.7
#' W <- discriminating_simililarity(X, k=length(labels)/2, sigma, cg)
#'
#' @references
#' Local similarity and diversity preserving discriminant projection for face and
#' handwriting digits recognition
#'
#' @export
discriminating_simililarity <- function(X, k=length(labels)/2, sigma, cg, threshold=.01) {
  #Wknn <- graph_weights(X)

  Wall <- graph_weights(X, k=k, weight_mode="heat", neighbor_mode="knn",sigma=sigma)
  Ww <- cg
  Wb <- 1- (cg > threshold)

  Ww2 <- Wall * Ww
  Wb2 <- Wall * Wb

  wind <- which(Ww2 >0)
  bind <- which(Wb2 >0)

  hw <- heat_kernel(Wall[wind], sigma)
  hb <- heat_kernel(Wall[bind], sigma)

  Wall[wind] <- hw * (1+hw)
  Wall[bind] <- hb * (1-hb)
  Wall
}



