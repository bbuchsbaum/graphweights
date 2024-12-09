#' Construct a Sparse Matrix of Spatial Constraints for Data Blocks
#'
#' This function creates a sparse matrix of spatial constraints for a set of data blocks. The spatial constraints matrix is useful in applications like image segmentation, where spatial information is crucial for identifying different regions in the image.
#'
#' @section Details:
#' The function computes within-block and between-block constraints based on the provided coordinates, bandwidths, and other input parameters. It then balances the within-block and between-block constraints using a shrinkage factor, and normalizes the resulting matrix by the first eigenvalue.
#'
#' @param coords The spatial coordinates as a matrix with rows as objects and columns as dimensions; or as a list of matrices where each element of the list contains the coordinates for a block.
#' @param nblocks The number of coordinate blocks. Default is 1.
#' @param sigma_within The bandwidth of the within-block smoother. Default is 5.
#' @param sigma_between The bandwidth of the between-block smoother. Default is 1.
#' @param shrinkage_factor The amount of shrinkage towards the spatial block average. Default is 0.1.
#' @param nnk_within The maximum number of nearest neighbors for within-block smoother. Default is 27.
#' @param nnk_between The maximum number of nearest neighbors for between-block smoother. Default is 1.
#' @param weight_mode_within The within-block nearest neighbor weight mode ("heat" or "binary"). Default is "heat".
#' @param weight_mode_between The between-block nearest neighbor weight mode ("heat" or "binary"). Default is "binary".
#' @param variable_weights A vector of per-variable weights. Default is 1.
#' @param verbose A boolean indicating whether to print progress messages. Default is FALSE.
#'
#' @return A sparse matrix representing the spatial constraints for the provided data blocks.
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
                                variable_weights=1, verbose=FALSE) {

  assert_that(shrinkage_factor > 0 & shrinkage_factor <= 1)


  if (is.list(coords)) {
    assert_that(length(coords) == nblocks)
    coords <- lapply(coords, as.matrix)
    nvars <- sum(sapply(coords, nrow))
  } else {
    coords <- as.matrix(coords)
    nvars <- nrow(coords) * nblocks
  }

  if (length(variable_weights) == 1) {
    variable_weights <- rep(1, nvars)
  } else {
    assert_that(length(variable_weights) == nvars)
  }


  if (verbose) {
    message("spatial_contraints: computing spatial adjacency (within)")
  }

  if (nblocks == 1) {
    Sw <- neighborweights::spatial_adjacency(coords, sigma=sigma_within,
                                           weight_mode=weight_mode_within,
                                           nnk=nnk_within,  normalized=FALSE,
                                           stochastic = TRUE)
    return(Sw)
  }

  Swithin <- if (length(coords) == 1) {
    if (verbose) {
      message("spatial_contraints: replicating within block constraints")
    }

    Matrix::bdiag(replicate(nblocks, Sw, simplify=FALSE))
    #indices <- rep(1:nrow(coords), nblocks)
  } else {
    if (verbose) {
      message("spatial_contraints: computing within-block constraints")
    }
    Sws <- lapply(1:nblocks, function(i) {
      neighborweights::spatial_adjacency(coords[[i]], sigma=sigma_within,
                                               weight_mode=weight_mode_within,
                                               nnk=nnk_within,  normalized=FALSE,
                                               stochastic = TRUE)
    })

    Matrix::bdiag(Sws)
  }

  Sbfin <- if (length(coords) == 1) {
    if (verbose) {
      message("spatial_contraints: computing spatial adjacency (between)")
    }

    Sb <- neighborweights::spatial_adjacency(coords,sigma=sigma_between,
                                           weight_mode=weight_mode_between,
                                           normalized=FALSE,
                                           nnk=nnk_between, stochastic=TRUE)


    Sbt <- as(Sb, "dgTMatrix")
    nvox <- nrow(coords)

    offsets <- cumsum(c(0, rep(nvox, nblocks-1)))

    out <- unlist(lapply(1:nblocks, function(i) {
      print(i)
      lapply(1:nblocks, function(j) {
      if (i == j) {
        NULL
      } else {
        cbind(Sbt@i + offsets[i] + 1, Sbt@j + offsets[j] +1, Sbt@x)
      }
      })
    }), recursive=FALSE)

    out <- out[!sapply(out, is.null)]
    out <- do.call(rbind, out)

    if (verbose) {
      message("spatial_contraints: building between matrix")
    }


    sparseMatrix(i=out[,1], j=out[,2], x=out[,3], dims=c(nvox*nblocks, nvox*nblocks))
  } else {
    allcds <- do.call(rbind, coords)
    index <- unlist(lapply(1:length(coords), function(i) rep(i, nrow(coords[[i]]))))
    b <- binary_label_matrix(index, index)
    Sb <- neighborweights::spatial_adjacency(allcds,sigma=sigma_between,
                                             weight_mode=weight_mode_between,
                                             normalized=FALSE,
                                             nnk=nnk_between)
    Sb[b > 0] <- 0
    Sb
  }

  if (verbose) {
    message("spatial_contraints: making doubly stochastic")
  }

  Sbfin <- neighborweights::make_doubly_stochastic(Sbfin)


  ## scale within matrix by variable weights
  if (any(variable_weights[1] != variable_weights)) {
    Wg <- Diagonal(x=sqrt(variable_weights))
    Swithin <- Wg %*% Swithin %*% Wg
  }

  if (verbose) {
    message("spatial_contraints: weighting by sfac")
  }


  ## compute ratio of within to between weights
  ## assumes all positive...
  rat <- sum(Swithin)/sum(Sbfin)
  sfac <- 1/shrinkage_factor * rat

  if (verbose) {
    message("spatial_contraints: balancing within and between by shrinkage_factor ", sfac)
  }

  ## balance within and between weights
  Stot <- (1-shrinkage_factor)*(1/rat)*Swithin + shrinkage_factor*Sbfin

  if (verbose) {
    message("spatial_contraints: normalizing by first eigenvalue")
  }

  S <- Stot/(RSpectra::eigs_sym(Stot, k=1, which="LA")$values[1])
  S

}

#' Construct Feature-Weighted Spatial Constraints for Data Blocks
#'
#' This function creates a sparse matrix of feature-weighted spatial constraints for a set of data blocks. The feature-weighted spatial constraints matrix is useful in applications like image segmentation and analysis, where both spatial and feature information are crucial for identifying different regions in the image.
#'
#' @section Details:
#' The function computes within-block and between-block constraints based on the provided coordinates, feature matrices, and other input parameters. It balances the within-block and between-block constraints using a shrinkage factor, and normalizes the resulting matrix by the first eigenvalue. The function also takes into account the weights of the variables in the provided feature matrices.
#'
#' @param coords The spatial coordinates as a matrix with rows as objects and columns as dimensions.
#' @param feature_mats A list of feature matrices, one for each data block.
#' @param sigma_within The bandwidth of the within-block smoother. Default is 5.
#' @param sigma_between The bandwidth of the between-block smoother. Default is 3.
#' @param wsigma_within The bandwidth of the within-block feature weights. Default is 0.73.
#' @param wsigma_between The bandwidth of the between-block feature weights. Default is 0.73.
#' @param alpha_within The scaling factor for within-block feature weights. Default is 0.5.
#' @param alpha_between The scaling factor for between-block feature weights. Default is 0.5.
#' @param shrinkage_factor The amount of shrinkage towards the spatial block average. Default is 0.1.
#' @param nnk_within The maximum number of nearest neighbors for within-block smoother. Default is 27.
#' @param nnk_between The maximum number of nearest neighbors for between-block smoother. Default is 27.
#' @param maxk_within The maximum number of nearest neighbors for within-block computation. Default is `nnk_within`.
#' @param maxk_between The maximum number of nearest neighbors for between-block computation. Default is `nnk_between`.
#' @param weight_mode_within The within-block nearest neighbor weight mode ("heat" or "binary"). Default is "heat".
#' @param weight_mode_between The between-block nearest neighbor weight mode ("heat" or "binary"). Default is "binary".
#' @param variable_weights A vector of per-variable weights. Default is a vector of ones with length equal to the product of the number of columns in the `coords` matrix and the length of `feature_mats`.
#' @param verbose A boolean indicating whether to print progress messages. Default is FALSE.
#'
#' @return A sparse matrix representing the feature-weighted spatial constraints for the provided data blocks.
#'
#' @examples
#' coords <- as.matrix(expand.grid(1:10, 1:10))
#' fmats <- replicate(20, matrix(rnorm(100*10), 10, 100), simplify=FALSE)
#' conmat <- feature_weighted_spatial_constraints(coords, fmats)
#'
#' conmat <- feature_weighted_spatial_constraints(coords, fmats, maxk_between=4, maxk_within=2,sigma_between=5, nnk_between=60)
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
                                                 variable_weights=rep(1, ncol(coords)*length(feature_mats)), verbose=FALSE) {

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
