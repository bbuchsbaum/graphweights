#' Local + Global KNN Adjacency
#'
#' Build an adjacency matrix that mixes \code{L} neighbors inside a local radius
#' \code{r} with \code{K} neighbors outside that radius. Far neighbors receive a
#' mild penalty so they can contribute without dominating.
#'
#' @param coord_mat Numeric matrix of coordinates (rows = points).
#' @param L Number of local neighbors (within \code{r}) to keep for each point.
#' @param K Number of far neighbors (outside \code{r}) to keep for each point.
#' @param r Radius defining the local ball.
#' @param weight_mode Weighting scheme name (e.g., "heat"); forwarded to internal helper get_neighbor_fun.
#' @param sigma Bandwidth for the heat/normalized kernels; default \code{r/2}.
#' @param far_penalty Either \code{"lambda"} (constant multiplier) or
#'   \code{"exp"} (decay with distance beyond \code{r}).
#' @param lambda Constant multiplier for far neighbors when
#'   \code{far_penalty = "lambda"}.
#' @param tau Scale of exponential decay when \code{far_penalty = "exp"}.
#' @param nnk_buffer Extra candidates requested from the NN search to ensure
#'   enough far neighbors are available.
#' @param include_diagonal Logical; keep self-loops.
#' @param symmetric Logical; if TRUE, symmetrize by averaging \code{A} and
#'   \code{t(A)}.
#' @param normalized Logical; if TRUE, row-normalize the matrix (stochastic).
#'
#' @return A sparse adjacency matrix mixing local and far neighbors.
#' @export
#'
#' @examples
#' set.seed(1)
#' coords <- matrix(runif(200), ncol = 2)
#' A <- local_global_adjacency(coords, L = 4, K = 3, r = 0.15,
#'                             weight_mode = "heat", lambda = 0.7)
#' Matrix::rowSums(A)[1:5]
local_global_adjacency <- function(coord_mat, L = 5, K = 5, r,
                                   weight_mode = c("heat", "binary"),
                                   sigma = r / 2,
                                   far_penalty = c("lambda", "exp"),
                                   lambda = 0.6, tau = r,
                                   nnk_buffer = 10,
                                   include_diagonal = FALSE,
                                   symmetric = TRUE, normalized = FALSE) {

  assertthat::assert_that(!missing(r))
  assertthat::assert_that(r > 0)
  assertthat::assert_that(L >= 0, K >= 0)

  weight_mode <- match.arg(weight_mode)
  far_penalty <- match.arg(far_penalty)

  n <- nrow(coord_mat)
  if (n < 2) {
    return(Matrix::sparseMatrix(i = integer(), j = integer(), x = numeric(), dims = c(n, n)))
  }

  k_req <- max(1, min(n - 1, L + K + nnk_buffer))

  nn <- Rnanoflann::nn(data = coord_mat, points = coord_mat, k = k_req)
  w_fun <- get_neighbor_fun(weight_mode, len = ncol(coord_mat), sigma = sigma)

  keep_trip <- lapply(seq_len(n), function(i) {
    idx <- nn$indices[i, ]
    d   <- nn$distances[i, ]

    ok <- !is.na(idx) & idx > 0 & idx != i & !is.na(d)
    idx <- idx[ok]
    d   <- d[ok]

    loc_mask <- d <= r
    loc_idx  <- idx[loc_mask]
    loc_d    <- d[loc_mask]
    if (length(loc_idx) > L) {
      loc_idx <- loc_idx[seq_len(L)]
      loc_d   <- loc_d[seq_len(L)]
    }

    far_idx  <- idx[!loc_mask]
    far_d    <- d[!loc_mask]
    if (length(far_idx) > K) {
      far_idx <- far_idx[seq_len(K)]
      far_d   <- far_d[seq_len(K)]
    }

    if (!include_diagonal && length(loc_idx) + length(far_idx) == 0) {
      return(NULL)
    }

    base_w <- w_fun(c(loc_d, far_d))
    far_adj <- if (far_penalty == "lambda") {
      rep(lambda, length(far_d))
    } else {
      exp(-(far_d - r) / tau)
    }

    w <- c(
      base_w[seq_along(loc_d)],
      base_w[seq_along(far_d) + length(loc_d)] * far_adj
    )

    data.frame(
      i = i,
      j = c(loc_idx, far_idx),
      x = w,
      stringsAsFactors = FALSE
    )
  })

  trip <- do.call(rbind, keep_trip)
  if (is.null(trip)) {
    return(Matrix::sparseMatrix(i = integer(), j = integer(), x = numeric(), dims = c(n, n)))
  }

  A <- Matrix::sparseMatrix(i = trip$i, j = trip$j, x = trip$x, dims = c(n, n))

  if (symmetric) {
    A <- (A + Matrix::t(A)) / 2
  }

  if (normalized) {
    rs <- Matrix::rowSums(A)
    invd <- 1 / pmax(rs, .Machine$double.eps)
    A <- Matrix::Diagonal(x = invd) %*% A
  }

  A
}
