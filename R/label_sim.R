
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
binary_label_matrix <- function(a, b=NULL, type=c("s", "d")) {
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
}


#' Create a Label Adjacency Matrix
#'
#' Constructs a label adjacency matrix based on two sets of labels `a` and `b`, creating edges depending on the `type` parameter and the similarity function `simfun`.
#'
#' @param a A vector of labels for the first set of data points.
#' @param b A vector of labels for the second set of data points.
#' @param type A character specifying the type of adjacency matrix to create, either "s" for same labels or "d" for different labels (default: "s").
#' @param simfun A function to determine the similarity between the labels (default: NULL). If NULL, the function will use "==" for type "s" and "!=" for type "d".
#' @param dim1 The dimension of the first set of data points (default: length(a)).
#' @param dim2 The dimension of the second set of data points (default: length(b)).
#'
#' @return A label adjacency matrix with edges determined by the label relationships between `a` and `b` and the similarity function `simfun`.
#'
#' @export
label_matrix2 <- function(a, b, type=c("s", "d"), simfun=NULL, dim1=length(a),
                          dim2=length(b)) {
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
    #ret[which(is.na(ret), arr.ind=TRUE)] <- 0
    ret

  }
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



#' Compute a Similarity or Distance Matrix for Categorical Labels
#'
#' Calculates a similarity or distance matrix between two categorical label vectors `a` and `b`.
#'
#' @param a The first categorical label vector.
#' @param b The second categorical label vector.
#' @param type The type of matrix to construct, either a similarity ('s') or distance ('d') matrix.
#' @param return_matrix A logical flag indicating whether to return the result as a sparse matrix (default: TRUE) or a triplet.
#' @param simfun An optional user-provided similarity function. If not provided, a default similarity function will be used.
#' @param dim1 The length of the first dimension.
#' @param dim2 The length of the second dimension.
#'
#' @return A similarity or distance matrix between the categorical label vectors `a` and `b`, either as a sparse matrix or a triplet.
#'
#' @importFrom Matrix sparseMatrix
#' @export
label_matrix <- function(a, b, type=c("s", "d"), return_matrix=TRUE, simfun=NULL , dim1=length(a), dim2=length(b)) {
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

  if (return_matrix) {
    sparseMatrix(i=out[,1], j=out[,2], x=out[,3], dims=c(dim1, dim2))
  } else {
    out
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


