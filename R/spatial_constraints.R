#' Construct a Sparse Matrix of Spatial Constraints for Data Blocks
#'
#' @description
#' Creates a sparse matrix of spatial constraints that captures the spatial relationships
#' between data points within and between blocks. This is particularly useful in applications
#' such as image segmentation, spatial clustering, and spatially-aware dimensionality reduction.
#'
#' @details
#' The function constructs a spatial constraint matrix by:
#' 1. Computing within-block spatial relationships using a Gaussian kernel
#' 2. Computing between-block spatial relationships (if nblocks > 1)
#' 3. Balancing within and between-block relationships using shrinkage
#' 4. Normalizing the final matrix by its largest eigenvalue
#'
#' The spatial relationships are computed based on:
#' * Within-block: Closer points have stronger connections
#' * Between-block: Connections between blocks are controlled by sigma_between
#' * Shrinkage: Controls the balance between within and between-block relationships
#'
#' Common use cases include:
#' * Image segmentation: Group pixels based on spatial proximity
#' * fMRI analysis: Model spatial relationships between brain regions
#' * Spatial clustering: Incorporate spatial constraints in clustering
#' * Dimensionality reduction: Preserve spatial structure in low-dimensional embeddings
#'
#' @param coords A matrix or list of matrices containing spatial coordinates:
#'   * Matrix: rows are objects, columns are spatial dimensions (e.g., x, y, z)
#'   * List: each element is a matrix of coordinates for one block
#' @param nblocks Integer specifying number of data blocks (default: 1)
#'   * 1: Single block analysis
#'   * >1: Multi-block analysis with between-block constraints
#' @param sigma_within Numeric bandwidth for within-block spatial kernel (default: 5)
#'   * Larger values: Smoother spatial relationships
#'   * Smaller values: More localized relationships
#' @param sigma_between Numeric bandwidth for between-block spatial kernel (default: 1)
#'   * Controls strength of connections between blocks
#'   * Smaller values create more separation between blocks
#' @param shrinkage_factor Numeric between 0 and 1 (default: 0.1)
#'   * Controls balance between within/between block constraints
#'   * 0: Only within-block constraints
#'   * 1: Maximum between-block influence
#' @param nnk_within Integer maximum nearest neighbors within blocks (default: 27)
#'   * Limits computational complexity
#'   * Should be large enough to ensure connectivity
#' @param nnk_between Integer maximum nearest neighbors between blocks (default: 1)
#'   * Controls sparsity of between-block connections
#'   * Smaller values create sparser connections
#' @param weight_mode_within Character specifying within-block weight type:
#'   * "heat": Gaussian kernel weights (default)
#'   * "binary": Binary weights for thresholded connections
#' @param weight_mode_between Character specifying between-block weight type:
#'   * "heat": Gaussian kernel weights
#'   * "binary": Binary weights (default)
#' @param variable_weights Numeric vector of weights per variable (default: 1)
#'   * Length must match total number of variables
#'   * Allows differential weighting of variables
#' @param verbose Logical for progress messages (default: FALSE)
#'
#' @return A sparse Matrix object representing spatial constraints where:
#'   * Positive values indicate spatial relationships
#'   * Larger values indicate stronger relationships
#'   * Matrix is normalized by its largest eigenvalue
#'
#' @examples
#' # Example 1: Simple 2D grid
#' coords <- expand.grid(x = 1:10, y = 1:10)
#' S <- spatial_constraints(as.matrix(coords), 
#'                         sigma_within = 2,
#'                         nnk_within = 8)
#'
#' # Example 2: Multi-block analysis
#' block1 <- expand.grid(x = 1:5, y = 1:5)
#' block2 <- expand.grid(x = 6:10, y = 6:10)
#' coords <- list(as.matrix(block1), as.matrix(block2))
#' S <- spatial_constraints(coords,
#'                         nblocks = 2,
#'                         sigma_within = 2,
#'                         sigma_between = 1,
#'                         shrinkage_factor = 0.2)
#'
#' # Example 3: 3D neuroimaging-like data
#' coords3d <- expand.grid(x = 1:5, y = 1:5, z = 1:5)
#' S <- spatial_constraints(as.matrix(coords3d),
#'                         sigma_within = 3,
#'                         nnk_within = 26,
#'                         weight_mode_within = "heat")
#'
#' @seealso
#' * \code{\link{feature_weighted_spatial_constraints}} for feature-aware constraints
#' * \code{\link{neighbor_graph}} for graph construction
#'
#' @importFrom assertthat assert_that
#' @importFrom Matrix sparseMatrix Diagonal bdiag
#' @importFrom RSpectra eigs_sym
#'
#' @export
spatial_constraints <- function(coords, nblocks=1,
                              sigma_within=5,
                              sigma_between=1,
                              shrinkage_factor=.1,
                              nnk_within=27,
                              nnk_between=1,
                              weight_mode_within="heat",
                              weight_mode_between="binary",
                              variable_weights=1, 
                              verbose=FALSE) {

  # Input validation
  assert_that(is.numeric(sigma_within) && sigma_within > 0,
              msg="sigma_within must be a positive number")
  assert_that(is.numeric(sigma_between) && sigma_between > 0,
              msg="sigma_between must be a positive number")
  assert_that(is.numeric(shrinkage_factor) && shrinkage_factor > 0 && shrinkage_factor <= 1,
              msg="shrinkage_factor must be between 0 and 1")
  assert_that(is.numeric(nnk_within) && nnk_within > 0 && nnk_within == floor(nnk_within),
              msg="nnk_within must be a positive integer")
  assert_that(is.numeric(nnk_between) && nnk_between > 0 && nnk_between == floor(nnk_between),
              msg="nnk_between must be a positive integer")
  assert_that(is.numeric(nblocks) && nblocks > 0 && nblocks == floor(nblocks),
              msg="nblocks must be a positive integer")
  assert_that(weight_mode_within %in% c("heat", "binary"),
              msg="weight_mode_within must be either 'heat' or 'binary'")
  assert_that(weight_mode_between %in% c("heat", "binary"),
              msg="weight_mode_between must be either 'heat' or 'binary'")
  assert_that(is.logical(verbose),
              msg="verbose must be a logical value")

  # Coords validation
  if (is.list(coords)) {
    assert_that(length(coords) == nblocks,
                msg="Length of coords list must match nblocks")
    # Check each matrix in the list
    coords <- lapply(coords, function(mat) {
      assert_that(is.matrix(mat) || is.data.frame(mat),
                 msg="Each element in coords must be a matrix or data frame")
      assert_that(all(sapply(mat, is.numeric)),
                 msg="All coordinates must be numeric")
      assert_that(nrow(mat) > 0,
                 msg="Empty coordinate matrices are not allowed")
      as.matrix(mat)
    })
    # Check consistent dimensionality
    dims <- unique(sapply(coords, ncol))
    assert_that(length(dims) == 1,
               msg="All coordinate matrices must have the same number of dimensions")
    nvars <- sum(sapply(coords, nrow))
  } else {
    assert_that(is.matrix(coords) || is.data.frame(coords),
               msg="coords must be a matrix, data frame, or list of matrices")
    assert_that(all(sapply(coords, is.numeric)),
               msg="All coordinates must be numeric")
    assert_that(nrow(coords) > 0,
               msg="Empty coordinate matrix is not allowed")
    coords <- as.matrix(coords)
    nvars <- nrow(coords) * nblocks
  }

  # Variable weights validation
  if (length(variable_weights) == 1) {
    assert_that(is.numeric(variable_weights) && variable_weights > 0,
               msg="variable_weights must be positive")
    variable_weights <- rep(variable_weights, nvars)
  } else {
    assert_that(length(variable_weights) == nvars,
               msg="Length of variable_weights must match total number of variables")
    assert_that(all(is.numeric(variable_weights) & variable_weights > 0),
               msg="All variable weights must be positive numbers")
  }

  # Check if nnk_within is smaller than the number of points in any block
  if (is.list(coords)) {
    min_points <- min(sapply(coords, nrow))
    assert_that(nnk_within < min_points,
               msg="nnk_within must be smaller than the number of points in each block")
  } else {
    assert_that(nnk_within < nrow(coords),
               msg="nnk_within must be smaller than the number of points")
  }

  if (verbose) {
    message("spatial_constraints: computing spatial adjacency (within)")
  }

  # Single block case
  if (nblocks == 1) {
    tryCatch({
      Sw <- neighborweights::spatial_adjacency(coords, sigma=sigma_within,
                                             weight_mode=weight_mode_within,
                                             nnk=nnk_within, normalized=FALSE,
                                             stochastic=TRUE)
      return(Sw)
    }, error = function(e) {
      stop("Error computing spatial adjacency: ", e$message)
    })
  }

  # Multi-block case
  Swithin <- if (length(coords) == 1) {
    if (verbose) {
      message("spatial_constraints: replicating within block constraints")
    }
    tryCatch({
      Matrix::bdiag(replicate(nblocks, Sw, simplify=FALSE))
    }, error = function(e) {
      stop("Error creating block diagonal matrix: ", e$message)
    })
  } else {
    if (verbose) {
      message("spatial_constraints: computing within-block constraints")
    }
    tryCatch({
      Sws <- lapply(1:nblocks, function(i) {
        neighborweights::spatial_adjacency(coords[[i]], sigma=sigma_within,
                                         weight_mode=weight_mode_within,
                                         nnk=nnk_within, normalized=FALSE,
                                         stochastic=TRUE)
      })
      Matrix::bdiag(Sws)
    }, error = function(e) {
      stop("Error computing within-block constraints: ", e$message)
    })
  }

  # Between-block computation
  Sbfin <- if (length(coords) == 1) {
    if (verbose) {
      message("spatial_constraints: computing spatial adjacency (between)")
    }

    tryCatch({
      Sb <- neighborweights::spatial_adjacency(coords, sigma=sigma_between,
                                             weight_mode=weight_mode_between,
                                             normalized=FALSE,
                                             nnk=nnk_between, stochastic=TRUE)

      Sbt <- as(Sb, "dgTMatrix")
      nvox <- nrow(coords)
      offsets <- cumsum(c(0, rep(nvox, nblocks-1)))

      out <- unlist(lapply(1:nblocks, function(i) {
        if (verbose) message(sprintf("Processing block %d of %d", i, nblocks))
        lapply(1:nblocks, function(j) {
          if (i == j) {
            NULL
          } else {
            cbind(Sbt@i + offsets[i] + 1, Sbt@j + offsets[j] + 1, Sbt@x)
          }
        })
      }), recursive=FALSE)

      out <- out[!sapply(out, is.null)]
      if (length(out) == 0) stop("No valid between-block connections found")
      out <- do.call(rbind, out)

      if (verbose) {
        message("spatial_constraints: building between matrix")
      }

      sparseMatrix(i=out[,1], j=out[,2], x=out[,3], 
                  dims=c(nvox*nblocks, nvox*nblocks))
    }, error = function(e) {
      stop("Error computing between-block constraints: ", e$message)
    })
  } else {
    tryCatch({
      allcds <- do.call(rbind, coords)
      index <- unlist(lapply(1:length(coords), function(i) rep(i, nrow(coords[[i]]))))
      b <- binary_label_matrix(index, index)
      Sb <- neighborweights::spatial_adjacency(allcds, sigma=sigma_between,
                                             weight_mode=weight_mode_between,
                                             normalized=FALSE,
                                             nnk=nnk_between)
      Sb[b > 0] <- 0
      Sb
    }, error = function(e) {
      stop("Error computing between-block constraints for multiple coordinate sets: ", e$message)
    })
  }

  if (verbose) {
    message("spatial_constraints: making doubly stochastic")
  }

  # Final matrix computations with error handling
  tryCatch({
    Sbfin <- neighborweights::make_doubly_stochastic(Sbfin)

    # Scale within matrix by variable weights
    if (any(variable_weights != 1)) {
      Wg <- Diagonal(x=sqrt(variable_weights))
      Swithin <- Wg %*% Swithin %*% Wg
    }

    if (verbose) {
      message("spatial_constraints: weighting by shrinkage factor")
    }

    # Compute ratio of within to between weights
    rat <- sum(Swithin)/sum(Sbfin)
    if (!is.finite(rat)) stop("Invalid ratio computed: weights sum to 0 or infinity")
    
    sfac <- 1/shrinkage_factor * rat

    if (verbose) {
      message(sprintf("spatial_constraints: balancing within and between by shrinkage_factor %g", sfac))
    }

    # Balance within and between weights
    Stot <- (1-shrinkage_factor)*(1/rat)*Swithin + shrinkage_factor*Sbfin

    if (verbose) {
      message("spatial_constraints: normalizing by first eigenvalue")
    }

    # Compute largest eigenvalue and normalize
    eig <- RSpectra::eigs_sym(Stot, k=1, which="LA")
    if (!is.finite(eig$values[1]) || eig$values[1] <= 0) {
      stop("Invalid eigenvalue computed: matrix may be singular or indefinite")
    }
    
    S <- Stot/eig$values[1]
    
    if (!all(is.finite(S@x))) {
      stop("Non-finite values in final matrix")
    }
    
    S
  }, error = function(e) {
    stop("Error in final matrix computations: ", e$message)
  })
}

#' Construct Feature-Weighted Spatial Constraints for Data Blocks
#'
#' @description
#' Creates a sparse matrix of spatial constraints that incorporates both spatial proximity
#' and feature similarity between data points. This is particularly useful when both
#' spatial location and feature characteristics are important for defining relationships
#' between data points.
#'
#' @details
#' The function combines spatial and feature information by:
#' 1. Computing spatial relationships using coordinate information
#' 2. Computing feature similarities using provided feature matrices
#' 3. Combining spatial and feature information using scaling parameters
#' 4. Balancing within and between-block relationships
#'
#' Feature weighting is controlled by:
#' * wsigma parameters: Control the sensitivity to feature differences
#' * alpha parameters: Balance between spatial and feature information
#'
#' Common applications include:
#' * Image segmentation with multiple channels
#' * Multi-modal neuroimaging analysis
#' * Spatial-temporal data analysis
#' * Pattern recognition in spatial data
#'
#' @param coords Matrix of spatial coordinates:
#'   * Rows: Data points
#'   * Columns: Spatial dimensions (e.g., x, y, z)
#' @param feature_mats List of feature matrices:
#'   * Each matrix contains features for one block
#'   * Rows correspond to spatial coordinates
#' @param sigma_within Numeric spatial bandwidth within blocks (default: 5)
#'   * Controls spatial smoothing within blocks
#' @param sigma_between Numeric spatial bandwidth between blocks (default: 3)
#'   * Controls spatial connections between blocks
#' @param wsigma_within Numeric feature bandwidth within blocks (default: 0.73)
#'   * Controls feature similarity sensitivity within blocks
#' @param wsigma_between Numeric feature bandwidth between blocks (default: 0.73)
#'   * Controls feature similarity sensitivity between blocks
#' @param alpha_within Numeric between 0 and 1 (default: 0.5)
#'   * 0: Only spatial information within blocks
#'   * 1: Only feature information within blocks
#' @param alpha_between Numeric between 0 and 1 (default: 0.5)
#'   * 0: Only spatial information between blocks
#'   * 1: Only feature information between blocks
#' @param shrinkage_factor Numeric between 0 and 1 (default: 0.1)
#'   * Controls balance between within/between block constraints
#' @param nnk_within Integer maximum neighbors within blocks (default: 27)
#' @param nnk_between Integer maximum neighbors between blocks (default: 27)
#' @param maxk_within Integer maximum neighbors for feature computation (default: nnk_within)
#' @param maxk_between Integer maximum neighbors for feature computation (default: nnk_between)
#' @param weight_mode_within Character weight type within blocks (default: "heat")
#' @param weight_mode_between Character weight type between blocks (default: "binary")
#' @param variable_weights Numeric vector of weights per variable (default: ones)
#' @param verbose Logical for progress messages (default: FALSE)
#'
#' @return A sparse Matrix object representing feature-weighted spatial constraints where:
#'   * Values represent combined spatial-feature relationships
#'   * Matrix is normalized by its largest eigenvalue
#'
#' @examples
#' # Example 1: Image segmentation with RGB channels
#' coords <- expand.grid(x = 1:10, y = 1:10)
#' # Create random RGB values for each pixel
#' features <- matrix(runif(100 * 3), 100, 3)
#' S <- feature_weighted_spatial_constraints(
#'   as.matrix(coords),
#'   list(features),
#'   sigma_within = 2,
#'   wsigma_within = 0.5,
#'   alpha_within = 0.3
#' )
#'
#' # Example 2: Multi-modal analysis
#' coords <- expand.grid(x = 1:5, y = 1:5)
#' features1 <- matrix(rnorm(25 * 4), 25, 4)  # First modality
#' features2 <- matrix(rnorm(25 * 3), 25, 3)  # Second modality
#' S <- feature_weighted_spatial_constraints(
#'   as.matrix(coords),
#'   list(features1, features2),
#'   sigma_within = 2,
#'   sigma_between = 1,
#'   alpha_within = 0.4,
#'   alpha_between = 0.6
#' )
#'
#' @seealso
#' * \code{\link{spatial_constraints}} for purely spatial constraints
#' * \code{\link{neighbor_graph}} for graph construction
#'
#' @importFrom assertthat assert_that
#' @importFrom Matrix sparseMatrix Diagonal bdiag
#' @importFrom RSpectra eigs_sym
#'
#' @export
feature_weighted_spatial_constraints <- function(coords,
                                              feature_mats,
                                              sigma_within=5,
                                              sigma_between=3,
                                              wsigma_within=.73,
                                              wsigma_between=.73,
                                              alpha_within=.5,
                                              alpha_between=.5,
                                              shrinkage_factor=.1,
                                              nnk_within=27,
                                              nnk_between=27,
                                              maxk_within=nnk_within,
                                              maxk_between=nnk_between,
                                              weight_mode_within="heat",
                                              weight_mode_between="binary",
                                              variable_weights=rep(1, ncol(coords)*length(feature_mats)), 
                                              verbose=FALSE) {

  assert_that(shrinkage_factor > 0 & shrinkage_factor <= 1)

  coords <- as.matrix(coords)
  nvox <- nrow(coords)
  nblocks <- length(feature_mats)

  Swl <- furrr::future_map(seq_along(feature_mats), function(i) {
    sw <- neighborweights::weighted_spatial_adjacency(coords, t(feature_mats[[i]]),
                                                      alpha=alpha_within,
                                                      wsigma=wsigma_within,
                                                      sigma=sigma_within,
                                                      weight_mode=weight_mode_within,
                                                      nnk=nnk_within, normalized=TRUE, stochastic=FALSE)
    neighborweights::make_doubly_stochastic(sw)
  })


  Swithin <- Matrix::bdiag(Swl)

  cmb <- t(combn(1:length(feature_mats), 2))

  offsets <- cumsum(c(0, rep(nvox, nblocks-1)))


  bet <- do.call(rbind, furrr::future_map(1:nrow(cmb), function(i) {
    a <- cmb[i,1]
    b <- cmb[i,2]
    sm <- neighborweights::cross_weighted_spatial_adjacency(
      coords, coords, t(feature_mats[[a]]), t(feature_mats[[b]]),
      wsigma=wsigma_between, weight_mode=weight_mode_between,
      alpha=alpha_between,
      nnk=nnk_between, maxk=maxk_between,
      sigma=sigma_between, normalized=FALSE)
    sm <- neighborweights::make_doubly_stochastic(sm)

    sm_nc <- as (sm, "dgTMatrix")
    r1 <- cbind (i = sm_nc@i + 1 + offsets[a], j = sm_nc@j + 1 + offsets[b], x = sm_nc@x)
    r2 <- cbind (i = sm_nc@j + 1 + offsets[b], j = sm_nc@i + 1 + offsets[a], x = sm_nc@x)
    rbind(r1,r2)
  }))


  Sbfin <- sparseMatrix(i=bet[,1], j=bet[,2], x=bet[,3], dims=c(nvox*nblocks, nvox*nblocks))

  Sbfin <- neighborweights::make_doubly_stochastic(Sbfin)


  ## scale within matrix by variable weights
  if (any(variable_weights[1] != variable_weights)) {
    Wg <- Diagonal(x=sqrt(variable_weights))
    Swithin <- Wg %*% Swithin %*% Wg
  }

  ## compute ratio of within to between weights
  rat <- sum(Swithin)/sum(Sbfin)
  sfac <- 1/shrinkage_factor * rat

  ## balance within and between weights
  Stot <- (1-shrinkage_factor)*(1/rat)*Swithin + shrinkage_factor*Sbfin

  S <- Stot/(RSpectra::eigs_sym(Stot, k=1, which="LA")$values[1])
  S

}
