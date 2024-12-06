#' Compute Temporal Autocorrelation Matrix
#'
#' @description
#' Computes the temporal autocorrelation matrix of a multivariate time series using
#' a sliding window approach. The function can optionally compute the inverse
#' correlation matrix using shrinkage estimation.
#'
#' @details
#' The function performs the following steps:
#' 1. Computes the correlation matrix of the input data
#' 2. Optionally inverts the correlation matrix using shrinkage estimation
#' 3. Constructs a band sparse matrix where each band represents correlations
#'    at different time lags up to the specified window size
#'
#' The resulting matrix captures the temporal dependencies in the data at different
#' time lags, which can be useful for time series analysis, signal processing,
#' and temporal pattern recognition.
#'
#' @param X A numeric matrix where rows represent observations and columns represent
#'   time points. Each row should contain a separate time series.
#' @param window Integer specifying the window size for computing autocorrelation.
#'   Must be greater than 1 and less than the number of columns in X. Controls how
#'   many time lags to consider. Default is 3.
#' @param inverse Logical indicating whether to compute the inverse of the correlation
#'   matrix using shrinkage estimation. This can be useful when the correlation
#'   matrix is ill-conditioned. Default is FALSE.
#'
#' @return A sparse symmetric Matrix object representing the temporal autocorrelation.
#'   The matrix has dimensions ncol(X) x ncol(X), where element (i,j) represents
#'   the average correlation between time points i and j within the specified window.
#'
#' @examples
#' \dontrun{
#' # Generate example multivariate time series data
#' n_series <- 10
#' n_timepoints <- 20
#' X <- matrix(rnorm(n_series * n_timepoints), 
#'            nrow = n_series, 
#'            ncol = n_timepoints)
#'
#' # Compute temporal autocorrelation with window size 3
#' autocor <- temporal_autocor(X, window = 3)
#'
#' # Compute inverse temporal autocorrelation
#' inv_autocor <- temporal_autocor(X, window = 3, inverse = TRUE)
#' }
#'
#' @seealso
#' * \code{\link{temporal_adjacency}} for computing temporal adjacency matrices
#' * \code{\link{temporal_laplacian}} for computing temporal Laplacian matrices
#'
#' @importFrom assertthat assert_that
#' @importFrom Matrix bandSparse
#' @importFrom corpcor invcor.shrink
#' @export
temporal_autocor <- function(X, window=3, inverse=FALSE) {
  # Input validation
  if (!is.matrix(X)) {
    stop("X must be a numeric matrix")
  }
  if (!is.numeric(X)) {
    stop("X must contain numeric values")
  }
  if (any(is.na(X))) {
    stop("X contains missing values (NA)")
  }
  if (ncol(X) < 3) {
    stop("X must have at least 3 columns (time points)")
  }
  if (!is.numeric(window) || length(window) != 1) {
    stop("window must be a single numeric value")
  }
  if (window <= 1 || window >= ncol(X)) {
    stop(sprintf("window must be > 1 and < ncol(X) (%d)", ncol(X)))
  }
  
  # Compute correlation matrix with error handling
  tryCatch({
    cmat <- cor(X)
  }, error = function(e) {
    stop("Failed to compute correlation matrix. Check for constant columns or other correlation issues.")
  })
  
  # Handle inverse computation with error checking
  if (inverse) {
    tryCatch({
      cmat <- corpcor::invcor.shrink(cmat)
    }, error = function(e) {
      stop("Failed to compute inverse correlation matrix using shrinkage estimation.")
    })
  }

  # Compute temporal structure
  rowmat <- matrix(rep(1:nrow(cmat), ncol(cmat)), nrow(cmat), ncol(cmat))
  colmat <- matrix(rep(1:ncol(cmat), each=nrow(cmat)), nrow(cmat), ncol(cmat))
  delta <- abs(rowmat - colmat)

  # Compute band values with error checking
  cvals <- tryCatch({
    sapply(1:window, function(i) {
      ind <- which(delta == i)
      if (length(ind) == 0) return(0)
      mean(cmat[ind], na.rm = TRUE)
    })
  }, error = function(e) {
    stop("Error computing temporal bands")
  })

  # Create sparse matrix
  bmat <- matrix(unlist(cvals), nrow(cmat), length(cvals), byrow=TRUE)
  bLis <- as.data.frame(bmat)
  A <- bandSparse(nrow(cmat), k = 1:window, diag = bLis, symmetric=TRUE)
  A
}

#' Compute Temporal Adjacency Matrix
#'
#' @description
#' Constructs a temporal adjacency matrix for a time series using a sliding window
#' approach. The adjacency can be binary or weighted using a heat kernel, capturing
#' the temporal relationships between time points.
#'
#' @details
#' The function creates an adjacency matrix where connections are established
#' between time points within a specified window. The strength of connections
#' can be either:
#' * Binary (1 for connected points, 0 otherwise)
#' * Heat kernel weighted (exponentially decaying with temporal distance)
#'
#' The resulting matrix is symmetric and includes self-connections (diagonal
#' elements are 1). This structure is useful for temporal clustering, smoothing,
#' and pattern detection in time series data.
#'
#' @param time A numeric vector representing time points. Values should be ordered
#'   and typically represent equally spaced time intervals.
#' @param weight_mode Character string specifying the weighting scheme:
#'   * "heat": Gaussian heat kernel weights (exponential decay with distance)
#'   * "binary": Binary weights (1 for points within window, 0 otherwise)
#'   Default is "heat".
#' @param sigma Numeric value controlling the width of the Gaussian kernel when
#'   weight_mode = "heat". Larger values create broader temporal connections.
#'   Only used when weight_mode = "heat". Default is 1.
#' @param window Integer specifying the size of the temporal window. Determines
#'   how many adjacent time points are connected. Must be >= 1. Default is 2.
#'
#' @return A sparse symmetric Matrix object representing temporal adjacency.
#'   The matrix has dimensions length(time) x length(time), where element (i,j)
#'   represents the connection strength between time points i and j.
#'
#' @examples
#' \dontrun{
#' # Create evenly spaced time points
#' time <- seq(0, 10, by = 0.5)
#'
#' # Compute binary temporal adjacency
#' adj_binary <- temporal_adjacency(time, 
#'                                 weight_mode = "binary",
#'                                 window = 2)
#'
#' # Compute heat kernel weighted temporal adjacency
#' adj_heat <- temporal_adjacency(time,
#'                               weight_mode = "heat",
#'                               sigma = 1.5,
#'                               window = 3)
#' }
#'
#' @seealso
#' * \code{\link{temporal_autocor}} for computing temporal autocorrelation
#' * \code{\link{temporal_laplacian}} for computing temporal Laplacian matrices
#'
#' @importFrom Matrix sparseMatrix t diag
#' @export
temporal_adjacency <- function(time, weight_mode = c("heat", "binary"), sigma=1, window=2) {
  # Input validation
  if (!is.numeric(time)) {
    stop("time must be a numeric vector")
  }
  if (length(time) < 2) {
    stop("time must have at least 2 points")
  }
  if (any(is.na(time))) {
    stop("time contains missing values (NA)")
  }
  if (!is.numeric(sigma) || length(sigma) != 1 || sigma <= 0) {
    stop("sigma must be a positive numeric value")
  }
  if (!is.numeric(window) || length(window) != 1 || window < 1) {
    stop("window must be a positive integer")
  }
  if (window > length(time)) {
    stop("window cannot be larger than the length of time series")
  }
  if (!is.ordered(time)) {
    warning("time points are not strictly ordered")
  }
  
  weight_mode <- match.arg(weight_mode)
  len <- length(time)

  # Define weight function with error handling
  wfun <- if (weight_mode == "binary") {
    function(x) 1
  } else if (weight_mode == "heat") {
    function(x) {
      w <- heat_kernel(x, sigma)
      if (any(is.na(w))) {
        stop("Error computing heat kernel weights")
      }
      w
    }
  } else {
    stop("Invalid weight_mode")
  }

  # Compute temporal relationships with error handling
  f <- function(t) {
    if (length(t) >= 1) {
      tryCatch({
        h <- wfun(sqrt((t[1] - t)^2))
        cbind(t[1], t, h)
      }, error = function(e) {
        stop("Error computing temporal weights: ", e$message)
      })
    } else {
      NULL
    }
  }

  # Create sparse matrix with error handling
  tryCatch({
    m <- do.call(rbind, runner(time, f, k=window, lag=0, simplify=TRUE))
    if (is.null(m) || nrow(m) == 0) {
      stop("Failed to compute temporal relationships")
    }
    
    sm <- sparseMatrix(i=m[,1], j=m[,2], x=m[,3], dims=c(len, len))
    sm <- (sm + t(sm))
    diag(sm) <- 1
    sm
  }, error = function(e) {
    stop("Failed to create adjacency matrix: ", e$message)
  })
}

#' Compute Temporal Laplacian Matrix
#'
#' @description
#' Constructs a temporal Laplacian matrix from a time series, representing the
#' difference between the degree matrix and the adjacency matrix. This matrix
#' captures the temporal structure and can be used for various time series
#' analysis tasks.
#'
#' @details
#' The temporal Laplacian is computed as L = D - A, where:
#' * D is the diagonal degree matrix (row sums of adjacency matrix)
#' * A is the temporal adjacency matrix
#'
#' The Laplacian matrix is useful for:
#' * Temporal clustering
#' * Dimensionality reduction
#' * Signal processing
#' * Pattern detection in time series
#'
#' The structure of the Laplacian depends on the chosen weight mode and parameters
#' used to construct the underlying adjacency matrix.
#'
#' @param time A numeric vector representing time points. Values should be ordered
#'   and typically represent equally spaced time intervals.
#' @param weight_mode Character string specifying the weighting scheme:
#'   * "heat": Gaussian heat kernel weights (exponential decay with distance)
#'   * "binary": Binary weights (1 for points within window, 0 otherwise)
#'   Default is "heat".
#' @param sigma Numeric value controlling the width of the Gaussian kernel when
#'   weight_mode = "heat". Larger values create broader temporal connections.
#'   Only used when weight_mode = "heat". Default is 1.
#' @param window Integer specifying the size of the temporal window. Determines
#'   how many adjacent time points are connected. Must be >= 1. Default is 2.
#'
#' @return A sparse symmetric Matrix object representing the temporal Laplacian.
#'   The matrix has dimensions length(time) x length(time). The Laplacian
#'   has special spectral properties useful for various analyses:
#'   * It is positive semi-definite
#'   * Its eigenvalues provide information about the temporal structure
#'   * The smallest non-zero eigenvalue indicates temporal connectivity
#'
#' @examples
#' \dontrun{
#' # Create evenly spaced time points
#' time <- seq(0, 10, by = 0.5)
#'
#' # Compute Laplacian with binary weights
#' lap_binary <- temporal_laplacian(time, 
#'                                 weight_mode = "binary",
#'                                 window = 2)
#'
#' # Compute Laplacian with heat kernel weights
#' lap_heat <- temporal_laplacian(time,
#'                               weight_mode = "heat",
#'                               sigma = 1.5,
#'                               window = 3)
#'
#' # Analyze eigenvalues for temporal structure
#' eigenvals <- eigen(lap_heat)$values
#' }
#'
#' @seealso
#' * \code{\link{temporal_adjacency}} for computing temporal adjacency matrices
#' * \code{\link{temporal_autocor}} for computing temporal autocorrelation
#'
#' @importFrom Matrix Diagonal rowSums
#' @export
temporal_laplacian <- function(time, weight_mode = c("heat", "binary"), sigma=1, window=2) {
  # Input validation
  if (!is.numeric(time)) {
    stop("time must be a numeric vector")
  }
  if (length(time) < 2) {
    stop("time must have at least 2 points")
  }
  if (!is.numeric(sigma) || length(sigma) != 1 || sigma <= 0) {
    stop("sigma must be a positive numeric value")
  }
  if (!is.numeric(window) || length(window) != 1 || window < 1) {
    stop("window must be a positive integer")
  }
  
  weight_mode <- match.arg(weight_mode)
  
  # Compute adjacency matrix with error handling
  adj <- tryCatch({
    temporal_adjacency(time, weight_mode, sigma, window)
  }, error = function(e) {
    stop("Failed to compute adjacency matrix: ", e$message)
  })
  
  # Compute Laplacian with error handling
  tryCatch({
    D <- Diagonal(x=rowSums(adj))
    if (any(is.na(D@x))) {
      stop("Invalid degree matrix (contains NA values)")
    }
    D - adj
  }, error = function(e) {
    stop("Failed to compute Laplacian matrix: ", e$message)
  })
}
