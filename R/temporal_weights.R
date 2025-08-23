#' Compute the temporal autocorrelation of a matrix
#'
#' This function computes the temporal autocorrelation of a given matrix using a specified window size and
#' optionally inverts the correlation matrix.
#'
#' @param X A numeric matrix for which to compute the temporal autocorrelation
#' @param window integer, the window size for computing the autocorrelation, must be between 1 and ncol(X) (default is 3)
#' @param inverse logical, whether to compute the inverse of the correlation matrix (default is FALSE)
#'
#' @return A sparse symmetric matrix representing the computed temporal autocorrelation
#'
#' @examples
#' # Create an example matrix
#' X <- matrix(rnorm(50), nrow = 10, ncol = 5)
#'
#' # Compute the temporal autocorrelation
#' result <- temporal_autocor(X, window = 2)
#'
#' @export
temporal_autocor <- function(X, window=3, inverse=FALSE) {
  assertthat::assert_that(window > 1 && window < ncol(X))
  cmat <- cor(X)
  if (inverse) {
    cmat <- corpcor::invcor.shrink(cmat)
  }

  rowmat <- matrix(rep(1:nrow(cmat), ncol(cmat)), nrow(cmat), ncol(cmat))
  colmat <- matrix(rep(1:ncol(cmat), each=nrow(cmat)), nrow(cmat), ncol(cmat))
  delta <- abs(rowmat - colmat)

  cvals <- sapply(1:window, function(i) {
    ind <- which(delta == i)
    mean(cmat[ind])
  })

  bmat <- matrix(unlist(cvals), nrow(cmat), length(cvals), byrow=TRUE)
  bLis <- as.data.frame(bmat)
  A <- bandSparse(nrow(cmat), k = 1:window, diagonals = bLis, symmetric=TRUE)
  A

  #fill cmat as distance matrix()
}


#' Compute the temporal adjacency matrix of a time series
#'
#' This function computes the temporal adjacency matrix of a given time series using a specified weight mode,
#' sigma, and window size.
#'
#' @param time A numeric vector representing a time series
#' @param weight_mode Character, the mode for computing weights, either "heat" or "binary" (default is "heat")
#' @param sigma Numeric, the sigma parameter for the heat kernel (default is 1)
#' @param window Integer, the window size for computing adjacency (default is 2)
#'
#' @return A sparse symmetric matrix representing the computed temporal adjacency
#'
#' @examples
#' # Create an example time series
#' time <- 1:10
#'
#' # Compute the temporal adjacency matrix using the heat weight mode
#' result <- temporal_adjacency(time, weight_mode = "heat", sigma = 1, window = 2)
#'
#' @export
temporal_adjacency <- function(time, weight_mode = c("heat", "binary"), sigma=1, window=2) {
  weight_mode <- match.arg(weight_mode)
  len <- length(time)

  wfun <- if (weight_mode == "binary") {
    function(x) 1
  } else if (weight_mode == "heat") {
    function(x) heat_kernel(x, sigma)
  } else {
    stop()
  }

  f <- function(t) {
    if (length(t) >= 1) {
      h <- wfun(sqrt((t[1] - t)^2))
      cbind(t[1], t, h)
    } else {
      NULL
    }
  }

  # Simple sliding window implementation since runner package is not available
  m <- NULL
  for (i in seq_along(time)) {
    end_idx <- min(i + window - 1, length(time))
    t_window <- time[i:end_idx]
    result <- f(t_window)
    if (!is.null(result)) {
      m <- rbind(m, result)
    }
  }

  # m <- do.call(rbind, rollApply(time, window=window, function(t) {
  #   if (length(t) >= 1) {
  #     h <- wfun(sqrt((t[1] - t)^2))
  #     cbind(t[1], t, h)
  #   } else {
  #     NULL
  #   }
  # }))

  sm <- sparseMatrix(i=m[,1], j=m[,2], x=m[,3], dims=c(len, len))
  sm <- (sm + t(sm))
  diag(sm) <- 1
  sm

}

#' Compute the temporal Laplacian matrix of a time series
#'
#' This function computes the temporal Laplacian matrix of a given time series using a specified weight mode,
#' sigma, and window size.
#'
#' @param time A numeric vector representing a time series
#' @param weight_mode Character, the mode for computing weights, either "heat" or "binary" (default is "heat")
#' @param sigma Numeric, the sigma parameter for the heat kernel (default is 1)
#' @param window Integer, the window size for computing adjacency (default is 2)
#'
#' @return A sparse symmetric matrix representing the computed temporal Laplacian
#'
#' @examples
#' # Create an example time series
#' time <- 1:10
#'
#' # Compute the temporal Laplacian matrix using the heat weight mode
#' result <- temporal_laplacian(time, weight_mode = "heat", sigma = 1, window = 2)
#'
#' @export
temporal_laplacian <- function(time, weight_mode = c("heat", "binary"), sigma=1, window=2) {
  weight_mode <- match.arg(weight_mode)
  adj <- temporal_adjacency(time, weight_mode, sigma, window)
  Diagonal(x=rowSums(adj))  - adj
}
