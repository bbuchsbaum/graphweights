

#' convolve_matrix
#'
#' @param X the data matrix
#' @param Kern the Kernel used to smooth the matrix
#' @importFrom Matrix Diagonal
#' @importFrom assertthat assert_that
#'
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

#' label_sim_matrix
#'
#' @param a the first categorical label vector
#' @param b the second categorical label vector
#' @param type construct a a similarity ('s') or distance ('d') matrix.
#' @param return_matrix return as a sparse matrix or a triplet.
#' @importFrom Matrix sparseMatrix
label_sim_matrix <- function(a, b, type=c("s", "d"), return_matrix=TRUE) {
  type <- match.arg(type)

  a.idx <- which(!is.na(a))
  b.idx <- which(!is.na(b))

  sfun <- function(x, y) {
    if (is.na(x) || is.na(y)) {
      0
    } else if (x == y) {
      1
    } else {
      0
    }
  }

  dfun <- function(x, y) {
    if (is.na(x) || is.na(y)) {
      0
    } else if (x == y) {
      0
    } else {
      1
    }
  }

  fun <- if (type == "s") sfun else dfun

  out <- lapply(a.idx, function(i) {
    ret <- lapply(b.idx, function(j) {
      if (fun(a[i], b[j])) {
        c(i,j,1)
      } else {
        NULL
      }
    })
    ret[sapply(ret, function(x) !is.null(x))]
  })

  out <- unlist(out, recursive=FALSE)
  out <- do.call(rbind, out)

  if (return_matrix) {
    sparseMatrix(i=out[,1], j=out[,2], x=out[,3], dims=c(length(a), length(b)))
  } else {
    out
  }
}
