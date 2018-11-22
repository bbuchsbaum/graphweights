

#' @export
binary_label_matrix <- function(a, b, type=c("s", "d"), dim1=length(a), dim2=length(b)) {
  assert_that(length(a) == length(b))

  levs <- unique(a)

  mlist <- list()
  for (i in seq_along(levs)) {
    lev <- levs[i]

    if (type == "s") {
      idx1 <- which(a == lev)
      idx2 <- which(b == lev)
    } else {
      idx1 <- which(a != lev)
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

  Reduce("+", mlist)
}


#' @export
label_matrix2 <- function(a, b, type=c("s", "d"), simfun=NULL, dim1=length(a), dim2=length(b)) {
  type <- match.arg(type)

  if (is.null(simfun) && type == "s") {
    simfun <- "=="
  } else if (is.null(simfun) && type == "d") {
    simfun <- "!="
  }


  if (type == "s") {
    out <- outer(a,b, simfun)
    ret <- Matrix(out, sparse = TRUE)
    if (any(is.na(ret))) {
      ret[which(is.na(ret), arr.ind=TRUE)] <- 0
    }

    ret <- ret * 1
    ret
  } else if (type == "d") {
    out <- outer(a,b, simfun)
    ret <- Matrix(out, sparse = TRUE)
    if (any(is.na(ret))) {
      ret[which(is.na(ret), arr.ind=TRUE)] <- 0
    }

    ret <- ret * 1
    #ret[which(is.na(ret), arr.ind=TRUE)] <- 0
    ret

  }
}



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



#' label_matrix
#'
#' @param a the first categorical label vector
#' @param b the second categorical label vector
#' @param type construct a a similarity ('s') or distance ('d') matrix.
#' @param return_matrix return as a sparse matrix or a triplet.
#' @param the length of the first dimension
#' @param the length of the second dimension
#' @importFrom Matrix sparseMatrix
#' @importFrom Matrix sparseMatrix
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


