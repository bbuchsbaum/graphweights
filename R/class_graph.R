#' @export
#' @param labels the class label vector
#' @importFrom Matrix sparseVector tcrossprod Matrix t
#'
#' @examples
#' data(iris)
#' labels <- iris[,5]
#' cg <- class_graph(labels)
#'
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


# local scatter
locality_preserving_scatter.class_graph <- function(x, ng, X) {
  dens <- node_density(ng)
  within_class_neighbors(x)
  ## need knn and class_graph
}


# equation (9,10)

#' @export
intraclass_scatter.class_graph <- function(x, X, q=1.5, mc.cores=1) {
  dens <-  intraclass_density(x, X, q, mc.cores)
  ipairs <- intraclass_pairs(x)
  ret <- lapply(1:nrow(ipairs), function(n) {
    i <- ipairs[n,1]
    j <- ipairs[n,2]
    xi <- X[i,]
    xj <- X[j,]
    dij <- sum( (xi - xj)^2)
    edij <- exp(-dij/dens[i])
    w <- 1/2 * edij * (1 + edij)
    cbind(i, j, w)
  })

  triplet <- do.call(rbind, ret)
  W <- sparseMatrix(i=triplet[,1], j=triplet[,2], x=triplet[,3])
  W <- (W + t(W))/2
}

# equation (9,10)

#' @export
interclass_scatter.class_graph <- function(x, X, q=1.5, mc.cores=1) {
  dens <-  interclass_density(x, X, q, mc.cores)
  ipairs <- interclass_pairs(x)

  ## need to go from i to j and j to i
  ret <- lapply(1:nrow(ipairs), function(n) {
    i <- ipairs[n,1]
    j <- ipairs[n,2]
    xi <- X[i,]
    xj <- X[j,]
    dij <- sum( (xi - xj)^2)
    edij <- exp(-dij/dens[i])
    w <- 1/2 * edij * (1 - edij)
    cbind(i, j, w)
  })

  triplet <- do.call(rbind, ret)
  W <- sparseMatrix(i=triplet[,1], j=triplet[,2], x=triplet[,3])
  W <- (W + t(W))/2
}


## Discriminative globality and locality preserving graph embedding for
## dimensionality reduction
## equations 11 (12)
#' @export
intraclass_density.class_graph <- function(x, X, q=1.5, mc.cores=1) {
  assert_that(nrow(X) == nvertices(x))
  labs <- x$labels

  #1/(n_c[i]^q) * sum((x_i - x_j)^2) for all x_i
  ret <- parallel::mclapply(x$levels, function(l) {
    idx <- which(labs == l)
    xc <- X[idx,]
    sapply(1:length(idx), function(i) {
      S <- sum(sweep(xc[-i,], 2, xc[i,], "-")^2)
      1/(length(idx)^q) * S
    })
  })

  dens <- numeric(nrow(X))

  for (i in 1:nclasses(x)) {
    dens[x$class_indices[[i]]] <- ret[[i]]
  }

  dens
}

## equation 18, 19
#' @export
interclass_density.class_graph <- function(x, X, q=2, mc.cores=1) {
  assert_that(nrow(X) == nvertices(x))
  labs <- x$labels

  vind <- 1:nvertices(x)
  sapply(1:nrow(X), function(i) {
    #nabes <- neighbors(x, i)[[1]]
    #non <- vind[vind %in% nabes]
    non <- vind[-match(nabes,vind)]
    #non <- which(!(x$data[i,,drop=FALSE] == 1))
    S <- sum(sweep(X[non,], 2, X[i,], "-")^2)
    1/(length(non)^q) * S
  })

}

#' @export
intraclass_pairs.class_graph <- function(x) {
  ## convert to triplet-style matrix?
  do.call(rbind, lapply(1:length(x$class_indices), function(i) {
    expand.grid(x=x$class_indices[[i]], y=x$class_indices[[i]])
  }))

}

#' @export
interclass_pairs.class_graph <- function(x) {
  do.call(rbind, lapply(1:length(x$class_indices), function(i) {
    expand.grid(x=x$class_indices[[i]], y=unlist(x$class_indices[-i]))
  }))

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

#' Compute similarity graph weighted by class structure
#'
#' @param X
#' @param k
#' @param sigma
#' @param cg
#'
#' @examples
#'
#' X <- matrix(rnorm(100*100), 100,100)
#' labels <- factor(rep(1:5, each=20))
#' cg <- class_graph(labels)
#' sigma <- .7
#'
#' @export
#'
## ref: Local similarity and diversity preserving discriminant projection for face and
## handwriting digits recognition
discriminating_simililarity <- function(X, k=length(labels)/2, sigma,cg, threshold=.01) {
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



