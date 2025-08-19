#' Create a Binary Label Adjacency Matrix
#'
#' Constructs a binary adjacency matrix based on two sets of labels `a` and `b`, creating edges when `a` and `b` have the same or different labels, depending on the `type` parameter.
#'
#' @param a A vector of labels for the first set of data points.
#' @param b A vector of labels for the second set of data points (default: NULL). If NULL, `b` will be set to `a`.
#' @param type A character specifying the type of adjacency matrix to create, either "s" for same labels or "d" for different labels (default: "s").
#'
#' @return A binary adjacency matrix with edges determined by the label relationships between `a` and `b`.
#'
#' @examples
#' data(iris)
#' a <- iris[,5]
#' bl <- binary_label_matrix(a, type="d")
#'
#' @export
#' @importFrom Matrix sparseVector tcrossprod
binary_label_matrix <- function(a, b = NULL, type = c("s", "d")) {
  type <- match.arg(type)
  if (is.null(b)) b <- a

  fa <- as.factor(a);  fb <- as.factor(b)
  na <- length(fa);    nb <- length(fb)

  if (type == "s") {
    # indices where levels coincide - need ALL pairs with same labels
    ia <- as.integer(fa)
    ib <- as.integer(fb)
    
    # Create all pairs (i,j) where label[i] == label[j]
    pairs <- expand.grid(i = 1:na, j = 1:nb)
    keep <- ia[pairs$i] == ib[pairs$j]
    
    sparseMatrix(i = pairs$i[keep],
                 j = pairs$j[keep],
                 x = 1L,
                 dims = c(na, nb))
  } else {
    # different labels - need ALL pairs with different labels
    ia <- as.integer(fa)
    ib <- as.integer(fb)
    
    # Create all pairs (i,j) where label[i] != label[j]
    pairs <- expand.grid(i = 1:na, j = 1:nb)
    keep <- ia[pairs$i] != ib[pairs$j]
    
    sparseMatrix(i = pairs$i[keep],
                 j = pairs$j[keep],
                 x = 1L,
                 dims = c(na, nb))
  }
}


#' Create a Label Adjacency Matrix
#'
#' Constructs a label adjacency matrix based on two sets of labels `a` and `b`, creating edges depending on the `type` parameter.
#'
#' @param a A vector of labels for the first set of data points.
#' @param b A vector of labels for the second set of data points.
#' @param type A character specifying the type of adjacency matrix to create, either "s" for same labels or "d" for different labels (default: "s").
#' @param dim1 The dimension of the first set of data points (default: length(a)).
#' @param dim2 The dimension of the second set of data points (default: length(b)).
#'
#' @return A label adjacency matrix with edges determined by the label relationships between `a` and `b` and the similarity function `simfun`.
#'
#' @export
label_matrix2 <- function(a, b, type = c("s", "d"), dim1 = length(a),
                          dim2 = length(b)) {

  type <- match.arg(type)
  fa <- as.factor(a);  fb <- as.factor(b)
  ia <- as.integer(fa); ib <- as.integer(fb)

  if (type == "s") {
    keep <- ia == ib
  } else {
    keep <- ia != ib
  }

  sparseMatrix(i = seq_len(dim1)[keep],
               j = seq_len(dim2)[keep],
               x = 1L,
               dims = c(dim1, dim2))
}



#' Convolve a Data Matrix with a Kernel Matrix
#'
#' Performs a convolution operation on a given data matrix `X` using a specified kernel matrix `Kern`.
#'
#' @param X A data matrix to be convolved.
#' @param Kern A kernel matrix used for the convolution.
#' @param normalize A logical flag indicating whether to normalize the kernel matrix before performing the convolution (default: FALSE).
#'
#' @return A matrix resulting from the convolution of the data matrix `X` with the kernel matrix `Kern`.
#'
#' @importFrom Matrix Diagonal
#' @importFrom assertthat assert_that
#' @export
convolve_matrix <- function(X, Kern, normalize = FALSE) {
  stopifnot(ncol(Kern) == nrow(Kern), ncol(X) == nrow(Kern))

  if (normalize) {
    rs <- rowSums(Kern)
    scale <- ifelse(rs > 0, 1 / sqrt(rs), 0)
    Kern <- Diagonal(x = scale) %*% Kern %*% Diagonal(x = scale)
  }
  X %*% Kern
}



#' Compute a Similarity or Distance Matrix for Categorical Labels
#'
#' Calculates a similarity or distance matrix between two categorical label vectors `a` and `b`.
#'
#' @param a The first categorical label vector.
#' @param b The second categorical label vector.
#' @param type The type of matrix to construct, either a similarity ('s') or distance ('d') matrix.
#' @param return_matrix A logical flag indicating whether to return the result as a sparse matrix (default: TRUE) or a triplet.
#' @param dim1 The length of the first dimension.
#' @param dim2 The length of the second dimension.
#'
#' @return A similarity or distance matrix between the categorical label vectors `a` and `b`, either as a sparse matrix or a triplet.
#'
#' @importFrom Matrix sparseMatrix
#' @export
label_matrix <- function(a, b, type = c("s", "d"),
                         return_matrix = TRUE,
                         dim1 = length(a),
                         dim2 = length(b)) {

  type <- match.arg(type)

  fa <- as.factor(a);  fb <- as.factor(b)
  ia <- as.integer(fa); ib <- as.integer(fb)

  if (type == "s") {
    keep <- ia == ib & !is.na(ia) & !is.na(ib)
  } else {                               # different
    keep <- ia != ib & !is.na(ia) & !is.na(ib)
  }

  if (!any(keep)) {
    if (return_matrix) {
      return(Matrix::Matrix(0, nrow = dim1, ncol = dim2, sparse = TRUE))
    } else {
      return(matrix(integer(0), ncol = 3))
    }
  }

  I <- which(keep)
  if (return_matrix) {
    sparseMatrix(i = I,
                 j = I,     # same length vectors ⇒ same indices
                 x = 1L,
                 dims = c(dim1, dim2))
  } else {
    cbind(I, I, 1L)
  }
}


#' Expand Similarity Between Labels Based on a Precomputed Similarity Matrix
#'
#' Expands the similarity between labels based on a precomputed similarity matrix, `sim_mat`, with either above-threshold or below-threshold values depending on the value of the `above` parameter.
#'
#' @param labels A vector of labels for which the similarities will be expanded.
#' @param sim_mat A precomputed similarity matrix containing similarities between the unique labels.
#' @param threshold A threshold value used to filter the expanded similarity values (default: 0).
#' @param above A boolean flag indicating whether to include the values above the threshold (default: TRUE) or below the threshold (FALSE).
#'
#' @return A sparse symmetric similarity matrix with the expanded similarity values.
#'
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

  # Ensure indices passed to C++ are integers and handle NAs by setting to 0 (or another indicator if needed)
  # C++ code currently skips negative indices, so 0 should work if C++ expects 1-based and subtracts 1.
  # Let's confirm C++ handles 0 correctly or adjust here.
  # C++ code does `indices[i] - 1`, so passing 0 will result in -1, which is skipped. This is okay.
  mind[is.na(mind)] <- 0
  mind <- as.integer(mind)

  out <- if (above) {
    expand_similarity_cpp(mind, sim_mat, threshold)
  } else {
    expand_similarity_below_cpp(mind, sim_mat, threshold)
  }

  # Check if the returned matrix is empty
  if (!is.matrix(out) || nrow(out) == 0) {
    warning("No similarities met the threshold criteria or C++ function failed.")
    # Return an empty sparse matrix with correct dimensions
    return(sparseMatrix(i={}, j={}, x={}, dims=c(length(labels), length(labels)), symmetric=TRUE))
  }

  # out is already a matrix, no need for do.call(rbind, out)
  # out <- do.call(rbind, out)

  # Ensure columns are numeric before passing to sparseMatrix
  i_col <- as.numeric(out[,1])
  j_col <- as.numeric(out[,2])
  x_col <- as.numeric(out[,3])
  
  # Check for NAs or non-finite values introduced somehow
  valid_idx <- is.finite(i_col) & is.finite(j_col) & is.finite(x_col)
  if (!all(valid_idx)) {
      warning("Non-finite values detected in matrix returned from C++; removing them.")
      i_col <- i_col[valid_idx]
      j_col <- j_col[valid_idx]
      x_col <- x_col[valid_idx]
  }
  
  sparseMatrix(i=i_col, j=j_col, x=x_col, 
               dims=c(length(labels), length(labels)), 
               symmetric=TRUE, 
               index1 = TRUE) # Assuming C++ returns 1-based indices

}


