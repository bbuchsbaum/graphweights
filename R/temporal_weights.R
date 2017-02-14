

#' temporal_adjacency
#'
#' @param time the time vector
#' @param weight_mode the type of weighting
#' @param sigma the bandwidth in units of time
#' @param window the window size in sampling units (assuming a regularily space grid)
#' @importFrom rowr rollApply
#' @export
temporal_adjacency <- function(time, weight_mode = c("heat", "binary"), sigma=1, window=2) {
  weight_mode <- match.arg(weight_mode)
  len <- length(time)

  wfun <- if (weight_mode == "binary") {
    function(x) 1
  } else {
    function(x) heat_kernel(x, sigma)
  }

  m <- do.call(rbind, rollApply(time, window=window, function(t) {
    if (length(t) > 1) {
      h <- wfun( (t[1] - t[-1])^2)
      cbind(t[1], t[-1], h)
    } else {
      NULL
    }
  }))


  sm <- sparseMatrix(i=m[,1], j=m[,2], x=m[,3], dims=c(len, len))
  (sm + t(sm))/2

}

#' temporal_laplacian
#' @param time
#' @param weight_mode
#' @param sigma
#' @param window
#' @export
temporal_laplacian <- function(time, weight_mode = c("heat", "binary"), sigma=1, window=2) {
  adj <- temporal_adjacency(time, weight_mode, sigma, window)
  Diagonal(x=rowSums(adj))  - adj
}
