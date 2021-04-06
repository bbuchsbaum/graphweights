#' temporal autocorrelation
#'
#' Compute a sparse temporal autocorrelation matrix for a temporal window
#'
#' @param X the data matrix, columns are variables, rows are time-points
#' @param window th number of time-points included in the auto-correlation band
#' @param inverse whether to return the inverse correlation matrix uisng `corpcor::invcor.shrink`
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

  bmat <- matrix(unlist(cvals), nrow(cmat), length(cvals) , byrow=TRUE)
  bLis <- as.data.frame(bmat)
  A <- bandSparse(nrow(cmat), k = 1:window, diag = bLis, symmetric=TRUE)
  A

  #fill cmat as distance matrix()
}



#' temporal_adjacency
#'
#' @param time the time vector
#' @param weight_mode the type of weighting
#' @param sigma the bandwidth in units of time
#' @param window the window size in sampling units (assuming a regularily space grid)
#' @export
#'
#' @import runner
#'
#' @examples
#'
#' time <- 1:100
#' tadj <- temporal_adjacency(time, weight_mode="heat", sigma=2, window=3)
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

  m  <- do.call(rbind, runner(time, f, k=window, lag=0, simplify=TRUE))

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

#' temporal_laplacian
#'
#' @inheritParams temporal_adjacency
#' @export
temporal_laplacian <- function(time, weight_mode = c("heat", "binary"), sigma=1, window=2) {
  weight_mode <- match.arg(weight_mode)
  adj <- temporal_adjacency(time, weight_mode, sigma, window)
  Diagonal(x=rowSums(adj))  - adj
}
