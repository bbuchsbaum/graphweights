#' Create a Label Similarity Matrix
#'
#' @description
#' Creates a similarity or dissimilarity matrix between categorical labels using the most
#' appropriate implementation based on the input characteristics and requirements.
#'
#' @details
#' This function provides a unified interface for creating label similarity matrices,
#' automatically selecting the most appropriate implementation based on:
#' * Input size and sparsity
#' * Memory constraints
#' * Custom similarity requirements
#'
#' Available implementations:
#' * Sparse: Memory-efficient for large, sparse label relationships
#' * Dense: Fast for smaller datasets with many relationships
#' * Custom: Flexible implementation supporting arbitrary similarity functions
#'
#' @param a Vector of labels for the first set of data points
#' @param b Optional vector of labels for the second set. If NULL, uses 'a'
#' @param type Type of relationship to compute:
#'   * "sparse": Memory-efficient implementation (formerly binary_label_matrix)
#'   * "dense": Fast implementation for smaller datasets (formerly label_matrix2)
#'   * "custom": Flexible implementation with custom similarity (formerly label_matrix)
#' @param similarity_fn Optional custom similarity function
#' @param threshold Threshold for similarity values (default: NULL)
#' @param mode Character, either "similarity" or "distance" (default: "similarity")
#'
#' @return A sparse Matrix object containing the computed similarities/distances
#'
#' @examples
#' # Basic sparse similarity matrix
#' labels <- factor(c("A", "B", "A", "C"))
#' sim_mat <- create_label_matrix(labels, type = "sparse")
#'
#' # Dense matrix with custom similarity
#' custom_sim <- function(x, y) as.numeric(x == y) * 2
#' sim_mat2 <- create_label_matrix(labels, type = "dense",
#'                                similarity_fn = custom_sim)
#'
#' @export
create_label_matrix <- function(a, b = NULL, 
                              type = c("sparse", "dense", "custom"),
                              similarity_fn = NULL,
                              threshold = NULL,
                              mode = c("similarity", "distance")) {
  
  type <- match.arg(type)
  mode <- match.arg(mode)
  
  # Convert mode to legacy type parameter
  legacy_type <- if (mode == "similarity") "s" else "d"
  
  # Select appropriate implementation
  result <- switch(type,
    "sparse" = .create_sparse_matrix(a, b, type = legacy_type),
    "dense" = .create_dense_matrix(a, b, type = legacy_type, 
                                 simfun = similarity_fn),
    "custom" = .create_custom_matrix(a, b, type = legacy_type, 
                                   simfun = similarity_fn)
  )
  
  # Apply threshold if specified
  if (!is.null(threshold)) {
    result@x[result@x < threshold] <- 0
    result <- drop0(result)
  }
  
  result
}

#' Create a Sparse Label Matrix (Internal)
#'
#' @keywords internal
.create_sparse_matrix <- function(a, b = NULL, type = c("s", "d")) {
  type <- match.arg(type)
  if (is.null(b)) {
    b <- a
  }

  if (is.factor(a)) {
    a <- as.character(a)
  }

  if (is.factor(b)) {
    b <- as.character(b)
  }

  assert_that(length(a) == length(b))

  levs <- unique(a)
  mlist <- list()

  for (i in seq_along(levs)) {
    lev <- levs[i]

    if (type == "s") {
      idx1 <- which(a == lev)
      idx2 <- which(b == lev)
    } else {
      idx1 <- which(a == lev)
      idx2 <- which(b != lev)
    }
    if (length(idx1) > 0 && length(idx2) > 0) {
      sv1=sparseVector(rep(1, length(idx1)), idx1, length(a))
      sv2=sparseVector(rep(1, length(idx2)), idx2, length(a))
      mlist[[i]] <- tcrossprod(sv1,sv2)
    }
  }

  mlist <- mlist[!sapply(mlist, is.null)]

  if (length(mlist) == 0) {
    stop("no overlapping levels in 'a' and 'b'")
  }

  ret <- Reduce("+", mlist)
  ret
}

#' Create a Dense Label Matrix (Internal)
#'
#' @keywords internal
.create_dense_matrix <- function(a, b, type = c("s", "d"), simfun = NULL) {
  type <- match.arg(type)

  if (is.null(simfun) && type == "s") {
    simfun <- "=="
  } else if (is.null(simfun) && type == "d") {
    simfun <- "!="
  }

  if (type == "s") {
    out <- outer(a,b, simfun)
    ret <- Matrix::Matrix(out, sparse = TRUE)

    if (any(is.na(ret))) {
      ret[which(is.na(ret), arr.ind=TRUE)] <- 0
    }

    ret <- ret * 1
    ret
  } else if (type == "d") {
    out <- outer(a,b, simfun)
    ret <- Matrix::Matrix(out, sparse = TRUE)
    if (any(is.na(ret))) {
      ret[which(is.na(ret), arr.ind=TRUE)] <- 0
    }

    ret <- ret * 1
    ret
  }
}

#' Create a Custom Label Matrix (Internal)
#'
#' @keywords internal
.create_custom_matrix <- function(a, b, type = c("s", "d"), simfun = NULL) {
  type <- match.arg(type)

  a.idx <- which(!is.na(a))
  b.idx <- which(!is.na(b))

  sfun <- if (!is.null(simfun)) {
    simfun
  } else {
    function(x, y) {
      if (is.na(x) || is.na(y)) {
        0
      } else if (x == y) {
        1
      } else {
        0
      }
    }
  }

  dfun <- if (!is.null(simfun)) {
    simfun
  } else {
    function(x, y) {
      if (is.na(x) || is.na(y)) {
        0
      } else if (x == y) {
        0
      } else {
        1
      }
    }
  }

  fun <- if (type == "s") sfun else dfun

  out <- lapply(a.idx, function(i) {
    ret <- lapply(b.idx, function(j) {
      v <- fun(a[i], b[j])
      if (v != 0) {
        c(i,j,v)
      }
    })

    ret[sapply(ret, function(x) !is.null(x))]
  })

  out <- unlist(out, recursive=FALSE)
  out <- do.call(rbind, out)

  sparseMatrix(i=out[,1], j=out[,2], x=out[,3], dims=c(length(a), length(b)))
}

#' Convolve a Data Matrix with a Kernel Matrix
#'
#' @description
#' Performs a convolution operation between a data matrix and a kernel matrix,
#' with optional kernel normalization. This function is useful for applying
#' smoothing or feature transformation operations in graph-based learning.
#'
#' @details
#' The convolution operation performs matrix multiplication between the input
#' data matrix and the kernel matrix: X * Kern. When normalize=TRUE, the kernel
#' is first normalized using the degree matrix, which can help prevent numerical
#' instabilities and ensure proper scaling.
#'
#' Performance characteristics:
#' * Time complexity: O(n*m*k) where n,m,k are matrix dimensions
#' * Space complexity: O(n*m) for the output matrix
#' * Most efficient when matrices are sparse
#'
#' Common use cases:
#' * Smoothing feature matrices using graph kernels
#' * Applying diffusion operations in semi-supervised learning
#' * Feature transformation in spectral methods
#'
#' @param X A data matrix to be convolved
#' @param Kern A kernel matrix used for the convolution
#' @param normalize Logical; whether to normalize the kernel matrix (default: FALSE)
#'
#' @return A matrix resulting from the convolution X * Kern
#'
#' @examples
#' # Create sample data
#' X <- matrix(rnorm(20), 4, 5)
#' Kern <- matrix(runif(25), 5, 5)
#' 
#' # Basic convolution
#' result <- convolve_matrix(X, Kern)
#' 
#' # Normalized convolution
#' result_norm <- convolve_matrix(X, Kern, normalize = TRUE)
#'
#' @seealso
#' * \code{\link{create_label_matrix}} for creating similarity matrices
#' * \code{\link{expand_label_similarity}} for expanding label relationships
#'
#' @importFrom Matrix Diagonal
#' @importFrom assertthat assert_that
#' @export
convolve_matrix <- function(X, Kern, normalize=FALSE) {
  assert_that(ncol(Kern) == nrow(Kern))
  assert_that(ncol(X) == nrow(Kern))
  
  if (normalize) {
    D1 <- Diagonal(x=1/sqrt(rowSums(Kern)))
    Kern <- D1 %*% Kern %*% D1
  }

  out <- X %*% Kern
  out
}

#' Expand Label Similarities Using a Pre-computed Similarity Matrix
#'
#' @description
#' Expands categorical label relationships using a pre-computed similarity matrix,
#' creating a sparse similarity graph that captures complex relationships between
#' labels. This function is particularly useful for semi-supervised learning and
#' transfer learning scenarios where you have prior knowledge about label relationships.
#'
#' @details
#' This function maps categorical labels to a similarity space defined by a pre-computed
#' similarity matrix. It's particularly useful when:
#' * Labels have hierarchical relationships (e.g., taxonomies)
#' * Labels have continuous similarity scores (e.g., word embeddings)
#' * You want to incorporate domain knowledge into label relationships
#'
#' Common use cases:
#' * Hierarchical classification with taxonomy relationships
#' * Text classification using word embedding similarities
#' * Transfer learning incorporating domain knowledge
#' * Semi-supervised learning with label propagation
#'
#' Performance characteristics:
#' * Time complexity: O(n^2) where n is the number of unique labels
#' * Space complexity: O(k) where k is the number of non-zero similarities
#' * Most efficient when the similarity matrix is sparse
#'
#' @param labels Vector of categorical labels to expand relationships for
#' @param sim_mat Pre-computed similarity matrix between unique labels.
#'   Row and column names must match the unique values in \code{labels}
#' @param threshold Numeric threshold for filtering similarities (default: 0)
#' @param above Logical; if TRUE keep similarities above threshold,
#'   if FALSE keep those below (default: TRUE)
#'
#' @return A sparse symmetric Matrix object where entry (i,j) represents the
#'   expanded similarity between labels[i] and labels[j]
#'
#' @examples
#' # Example with hierarchical categories
#' categories <- c("animal", "mammal", "cat", "dog", "plant", "tree")
#' n <- length(categories)
#' 
#' # Create a similarity matrix based on hierarchy
#' sim_mat <- matrix(0, n, n)
#' colnames(sim_mat) <- rownames(sim_mat) <- categories
#' 
#' # Define hierarchical relationships
#' sim_mat["mammal", "animal"] <- sim_mat["animal", "mammal"] <- 0.8
#' sim_mat["cat", "mammal"] <- sim_mat["mammal", "cat"] <- 0.9
#' sim_mat["dog", "mammal"] <- sim_mat["mammal", "dog"] <- 0.9
#' sim_mat["tree", "plant"] <- sim_mat["plant", "tree"] <- 0.9
#' 
#' # Sample data with these categories
#' labels <- c("cat", "dog", "tree", "plant", "animal")
#' 
#' # Expand similarities
#' expanded <- expand_label_similarity(labels, sim_mat, threshold = 0.5)
#'
#' @seealso
#' * \code{\link{create_label_matrix}} for creating basic label similarity matrices
#'
#' @importFrom Matrix Diagonal
#' @importFrom assertthat assert_that
#' @export
expand_label_similarity <- function(labels, sim_mat, threshold=0, above=TRUE) {
  cnames <- colnames(sim_mat)
  rnames <- rownames(sim_mat)
  assertthat::assert_that(!(is.null(cnames) && is.null(rnames)))

  if (!is.null(cnames) && !is.null(rnames)) {
    assertthat::assert_that(all.equal(cnames, rnames))
  }
  if (is.null(cnames)) {
    lnames <- rnames
  } else {
    lnames <- cnames
  }

  mind <- match(labels, lnames)

  if (all(is.na(mind))) {
    stop(paste("no matches between `labels` and similarity matrix entries"))
  }

  mind[is.na(mind)] <- 0

  out <- if (above) {
    expand_similarity_cpp(mind, sim_mat, threshold)
  } else {
    expand_similarity_below_cpp(mind, sim_mat, threshold)
  }

  if (length(out) == 0) {
    stop("similarity matching failed: no returned entries")
  }

  out <- do.call(rbind, out)
  sparseMatrix(i=out[,1], j=out[,2], x=out[,3], dims=c(length(labels), length(labels)), symmetric=TRUE)
}
