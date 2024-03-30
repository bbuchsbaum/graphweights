
# Algorithm 1: DGLPGE Method

# Step 1: Calculate the weight matrix W_plus of the global intra-class graph G_plus
W_plus <- compute_W_plus(X, labels)

# Step 1 continued: Calculate the intra-class globality preserving scatter S_W
S_W <- compute_S_W(X, W_plus)

# Step 2: Calculate the weight matrix W_minus of the global inter-class graph G_minus
W_minus <- compute_W_minus(X, labels)

# Step 2 continued: Calculate the inter-class globality preserving scatter S_B
S_B <- compute_S_B(X, W_minus)

# Step 3: Calculate the weight matrix W_k of the local graph G_k
W_k <- compute_W_k(X, k, labels)

# Step 3 continued: Calculate the locality preserving scatter S_L
S_L <- compute_S_L(X, W_k)

# Step 4: Solve Eq. (35) to form optimal projection P
P <- compute_optimal_projection(X, S_W, S_B, S_L, alpha, lambda)

# Now, you can use the optimal projection P to transform your data
Y <- X %*% P



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

compute_intraclass_scatter.plus <- function(x, data) {
  labels <- as.factor(x$labels)
  dist_data <- dist(data)
  dist_matrix <- as.matrix(dist_data)

  out <- Reduce("+", lapply(levels(labels), function(lev) {
    same_class_indices <- labels == lev
    dist_same_class <- dist_matrix[same_class_indices, same_class_indices]
    w_same_class <- 0.5 * (exp(-dist_same_class^2) + exp(-t(dist_same_class)^2))
    kronecker(Matrix(same_class_indices, sparse=FALSE), t(Matrix(same_class_indices, sparse=FALSE))) * w_same_class
  }))

  out
}

compute_W_plus <- function(class_graph, q = 1) {
  n <- nrow(class_graph$X)
  W_plus <- matrix(0, n, n)

  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      if (class_graph$labels[i] == class_graph$labels[j]) {
        t_i_plus <- sum((class_graph$X[i,] - class_graph$X[class_graph$class_indices[[class_graph$labels[i]]],])^2) / (length(class_graph$class_indices[[class_graph$labels[i]]])^q)
        t_j_plus <- sum((class_graph$X[j,] - class_graph$X[class_graph$class_indices[[class_graph$labels[j]]],])^2) / (length(class_graph$class_indices[[class_graph$labels[j]]])^q)

        w_ij_plus <- 0.5 * exp(-sum((class_graph$X[i,] - class_graph$X[j,])^2) / t_i_plus) * (1 + exp(-sum((class_graph$X[i,] - class_graph$X[j,])^2) / t_i_plus))
        w_ji_plus <- 0.5 * exp(-sum((class_graph$X[j,] - class_graph$X[i,])^2) / t_j_plus) * (1 + exp(-sum((class_graph$X[j,] - class_graph$X[i,])^2) / t_j_plus))

        W_plus[i, j] <- 0.5 * (w_ij_plus + w_ji_plus)
        W_plus[j, i] <- W_plus[i, j]
      }
    }
  }

  return(W_plus)
}

#' @export
compute_W_plus.class_graph <- function(class_graph, q = 1) {
  compute_W_plus(class_graph, q)
}

compute_W_minus <- function(class_graph, q = 1) {
  n <- nrow(class_graph$X)
  W_minus <- matrix(0, n, n)

  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      if (class_graph$labels[i] != class_graph$labels[j]) {
        t_i_minus <- sum((class_graph$X[i,] - class_graph$X[class_graph$class_indices[[class_graph$labels[i]]],])^2) / (length(class_graph$class_indices[[class_graph$labels[i]]])^q)
        t_j_minus <- sum((class_graph$X[j,] - class_graph$X[class_graph$class_indices[[class_graph$labels[j]]],])^2) / (length(class_graph$class_indices[[class_graph$labels[j]]])^q)

        w_ij_minus <- 0.5 * exp(-sum((class_graph$X[i,] - class_graph$X[j,])^2) / t_i_minus) * (1 - exp(-sum((class_graph$X[i,] - class_graph$X[j,])^2) / t_i_minus))
        w_ji_minus <- 0.5 * exp(-sum((class_graph$X[j,] - class_graph$X[i,])^2) / t_j_minus) * (1 - exp(-sum((class_graph$X[j,] - class_graph$X[i,])^2) / t_j_minus))

        W_minus[i, j] <- 0.5 * (w_ij_minus + w_ji_minus)
        W_minus[j, i] <- W_minus[i, j]
      }
    }
  }

  return(W_minus)
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

compute_Sw <- function(class_graph) {
  # Extract necessary components from the class_graph object
  W_plus <- class_graph$W
  X <- class_graph$X

  # Calculate the diagonal matrix D_plus
  D_plus <- Matrix::diag(Matrix::rowSums(W_plus))

  # Calculate the Laplacian matrix L_plus
  L_plus <- D_plus - W_plus

  # Compute the intraclass globality preserving scatter Sw
  Sw <- X %*% L_plus %*% t(X)

  # Return the computed Sw matrix
  return(Sw)
}

# Define the S3 method for class_graph objects
#' @export
compute_Sw.class_graph <- function(class_graph) {
  compute_Sw(class_graph)
}




# Function to compute Sb
compute_Sb <- function(class_graph) {
  W_minus <- compute_W_minus(class_graph)
  X <- class_graph$X

  # Calculate the diagonal matrix D_minus
  D_minus <- Matrix::diag(Matrix::rowSums(W_minus))

  # Calculate the Laplacian matrix L_minus
  L_minus <- D_minus - W_minus

  # Compute the inter-class globality preserving scatter Sb
  Sb <- X %*% L_minus %*% t(X)

  # Return the computed Sb matrix
  return(Sb)
}

# Define the S3 methods for class_graph objects
#' @export
compute_W_minus.class_graph <- function(class_graph) {
  compute_W_minus(class_graph)
}

#' @export
compute_Sb.class_graph <- function(class_graph) {
  compute_Sb(class_graph)
}

# Function to compute W_k
compute_Wk <- function(class_graph, k) {
  W_plus <- class_graph$W
  W_minus <- compute_W_minus(class_graph)

  # Compute the k-nearest neighbors for each point
  knn_indices <- classGraph::knn_index(class_graph, k)

  # Initialize the W_k matrix
  W_k <- Matrix(0, nrow = nrow(W_plus), ncol = ncol(W_plus))

  for (i in seq_along(class_graph$labels)) {
    for (j in knn_indices[[i]]) {
      if (class_graph$labels[i] == class_graph$labels[j]) {
        # Same class
        W_k[i, j] <- 0.5 * (W_plus[i, j] + W_plus[j, i])
      } else {
        # Different class
        W_k[i, j] <- -0.5 * (W_minus[i, j] + W_minus[j, i])
      }
    }
  }

  # Return the computed W_k matrix
  return(W_k)
}

# Function to compute Sl
compute_Sl <- function(class_graph, k) {
  W_k <- compute_Wk(class_graph, k)
  X <- class_graph$X

  # Calculate the diagonal matrix D_k
  D_k <- Matrix::diag(Matrix::rowSums(W_k))

  # Calculate the Laplacian matrix L_k
  L_k <- D_k - W_k

  # Compute the locality preserving scatter Sl
  Sl <- X %*% L_k %*% t(X)

  # Return the computed Sl matrix
  return(Sl)
}

# Define the S3 methods for class_graph objects
#' @export
compute_Wk.class_graph <- function(class_graph, k) {
  compute_Wk(class_graph, k)
}

#' @export
compute_Sl.class_graph <- function(class_graph, k) {
  compute_Sl(class_graph, k)
}

library(geigen)

solve_optimal_projection <- function(class_graph, alpha, k) {
  # Compute the matrices L_plus, L_minus, and L_k
  L_plus <- compute_L_plus(class_graph)
  L_minus <- compute_L_minus(class_graph)
  L_k <- compute_Lk(class_graph, k)

  # Compute the combined Laplacian matrix
  L_combined <- alpha * (L_minus - L_plus) + (1 - alpha) * L_k

  # Get the data matrix X
  X <- class_graph$X

  # Compute the matrix A
  A <- X %*% L_combined %*% t(X)

  # Solve the generalized eigenvalue problem
  eig_result <- geigen::geigen(A, X)

  # Sort the eigenvectors by their corresponding eigenvalues
  P <- eig_result$vectors[, order(eig_result$values)]

  # Return the optimal projection matrix P
  return(P)
}

#' @export
solve_optimal_projection.class_graph <- function(class_graph, alpha, k) {
  solve_optimal_projection(class_graph, alpha, k)
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

