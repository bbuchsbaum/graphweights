
# Design similarity kernel functions
# 
# A toolkit to build flexible design-similarity kernels for factorial designs,
# with support for nominal/ordinal factors, interactions, and mapping
# from cell space to effect-coded space.

#' @importFrom Matrix kronecker crossprod diag t Matrix
#' @importFrom stats contr.sum contr.helmert
NULL

#' Build a flexible design-similarity kernel
#' 
#' Constructs a positive semi-definite (PSD) kernel K over design regressors
#' that encodes factorial similarity: same level similarity, optional smoothness
#' across ordinal levels, and interaction structure via Kronecker composition.
#' Works either directly in cell space (one regressor per cell) or in
#' effect-coded space (main effects, interactions), by pulling K back
#' through user-specified contrast matrices.
#'
#' @details
#' The kernel is constructed using the following principles:
#' \itemize{
#'   \item For each factor i with L_i levels, define a per-factor kernel K_i (nominal or ordinal)
#'   \item For any term S (subset of factors), construct a term kernel by Kronecker product:
#'         K_S = kronecker(K_i, K_j, ...) for i,j in S and kronecker with J matrices for factors not in S
#'   \item Combine term kernels with nonnegative weights rho[S] and optionally add a small ridge
#'   \item If using effect coding, map the cell kernel to effect-space: K_effect = T' K_cell T
#' }
#'
#' @param factors A named list, one entry per factor, e.g.
#'   list(A=list(L=5,type="nominal"), B=list(L=5,type="ordinal", l=1.5)). For
#'   ordinal, supply length-scale l (>0). Supported types: "nominal", "ordinal",
#'   "circular" (wrap-around distances).
#' @param terms A list of character vectors; each character vector lists factor
#'   names that participate in that term. Examples: list("A", "B", c("A","B"))
#'   for main A, main B, A:B. If NULL, defaults to all singletons (main effects)
#'   and the full interaction.
#' @param rho A named numeric vector of nonnegative weights for each term (names
#'   like "A", "B", "A:B"). If NULL, defaults to 1 for each term. Use rho0 for
#'   the identity term (see below).
#' @param include_intercept Logical; if TRUE, adds rho0 * I to K (small ridge /
#'   identity term).
#' @param rho0 Nonnegative scalar weight for the identity term; default 1e-8.
#' @param basis Either "cell" (default) or "effect". If "effect", you must
#'   supply `contrasts`.
#' @param contrasts A named list of contrast matrices for each factor, with
#'   dimensions L_i x d_i. For example: list(A = contr.sum(5), B = contr.sum(5)).
#'   You can also pass orthonormal Helmert.
#' @param block_structure If basis="effect", a character vector naming the
#'   sequence of effect blocks to include, e.g., c("A","B","A:B"). If NULL,
#'   inferred from `terms` and `contrasts`.
#' @param normalize One of c("none","unit_trace","unit_fro","max_diag").
#' @param jitter Small diagonal added to ensure SPD; default 1e-8.
#' @return A list with elements:
#'   \item{K}{PSD kernel matrix in the requested basis (cell or effect)}
#'   \item{K_cell}{The cell-space kernel matrix (always returned)}
#'   \item{info}{A list containing levels, factor_names, term_names, basis, map
#'     (T matrix for effect basis), blocks}
#' @examples
#' # Simple 2x3 design with nominal factors
#' factors <- list(
#'   A = list(L=2, type="nominal"),
#'   B = list(L=3, type="nominal")
#' )
#' K1 <- design_kernel(factors)
#' print(dim(K1$K))  # 6x6 cell-space kernel
#' 
#' # Ordinal factor with smoothness
#' factors <- list(
#'   dose = list(L=3, type="ordinal", l=1.0),
#'   treat = list(L=2, type="nominal")
#' )
#' K2 <- design_kernel(factors)
#' print(dim(K2$K))  # 6x6 cell-space kernel
#' @export
design_kernel <- function(factors,
                          terms = NULL,
                          rho = NULL,
                          include_intercept = TRUE,
                          rho0 = 1e-8,
                          basis = c("cell","effect"),
                          contrasts = NULL,
                          block_structure = NULL,
                          normalize = c("none","unit_trace","unit_fro","max_diag"),
                          jitter = 1e-8) {
  
  # helper: %||%
  `%||%` <- function(a, b) if (is.null(a)) b else a
  
  basis <- match.arg(basis)
  normalize <- match.arg(normalize)

  # --- 1) per-factor kernels and utilities ---
  # Input validation
  if (length(factors) == 0) stop("`factors` must be a non-empty list.")
  fact_names <- names(factors)
  if (is.null(fact_names) || any(!nzchar(fact_names))) stop("`factors` must be a NAMED list.")
  if (length(unique(fact_names)) != length(fact_names)) stop("Factor names must be unique.")
  
  Ls <- sapply(factors, function(f) { 
    if (is.null(f$L)) stop("Each factor needs L (levels).") 
    if (!is.numeric(f$L) || f$L < 1) stop("Factor levels must be positive integers.")
    as.integer(f$L)
  })
  
  types <- sapply(factors, function(f) {
    type <- tolower(f$type %||% "nominal")
    if (!type %in% c("nominal", "ordinal", "circular")) {
      stop("Factor type must be 'nominal', 'ordinal', or 'circular'.")
    }
    type
  })
  
  lscl <- sapply(seq_along(factors), function(i) {
    if (types[i] %in% c("ordinal", "circular")) {
      l <- factors[[i]]$l %||% 1.0
      if (!is.numeric(l) || l <= 0) stop("Length scale l must be positive for ordinal/circular factors.")
      l
    } else {
      NA_real_
    }
  })

  # level index helpers
  make_levels <- function(L) 1:L

  # per-factor base kernel
  k_factor <- function(L, type="nominal", l=1.0) {
    idx <- make_levels(L)
    if (type == "nominal") {
      return(diag(L))
    } else if (type == "ordinal") {
      D <- outer(idx, idx, function(i,j) (i-j)^2)
      return(exp(- D / (2*l*l)))
    } else if (type == "circular") {
      # wrap-around distances on a ring of L
      D <- outer(idx, idx, function(i,j) {
        d <- abs(i-j)
        pmin(d, L - d)^2
      })
      return(exp(- D / (2*l*l)))
    } else {
      stop("Unknown factor type: ", type)
    }
  }

  J_list <- lapply(Ls, function(L) matrix(1, L, L))
  K_list <- mapply(function(L, type, l) k_factor(L, type, ifelse(is.na(l), 1.0, l)),
                   Ls, types, lscl, SIMPLIFY = FALSE)
  names(K_list) <- fact_names

  # --- 2) default terms if none specified ---
  if (is.null(terms)) {
    # singletons + full interaction
    terms <- c(as.list(fact_names), list(fact_names))
  }
  term_names <- vapply(terms, function(S) paste(S, collapse=":"), character(1))

  # --- 3) assemble K_cell by ANOVA composition ---
  # K_cell = sum_{S in terms} rho[S] * (⊗_{i in S} K_i ⊗ ⊗_{i not in S} J_i) + rho0 * I
  if (is.null(rho)) rho <- setNames(rep(1.0, length(terms)), term_names)
  if (any(rho < 0)) stop("All rho must be nonnegative.")
  if (is.null(names(rho))) names(rho) <- term_names

  kron_all <- function(mats) {
    out <- mats[[1]]
    if (length(mats) > 1) {
      for (i in 2:length(mats)) {
        # Use Matrix package kronecker for potential sparsity
        out <- Matrix::kronecker(out, mats[[i]])
      }
    }
    out
  }

  per_term_kron <- function(S) {
    mats <- vector("list", length(fact_names))
    for (i in seq_along(fact_names)) {
      nm <- fact_names[i]
      mats[[i]] <- if (nm %in% S) K_list[[nm]] else J_list[[i]]
    }
    kron_all(mats)
  }

  # size of cell space
  Qcell <- prod(Ls)
  K_cell <- matrix(0, Qcell, Qcell)
  for (k in seq_along(terms)) {
    S <- terms[[k]]
    K_cell <- K_cell + (rho[[ term_names[k] ]] %||% 1.0) * per_term_kron(S)
  }
  if (include_intercept && rho0 > 0) {
    K_cell <- K_cell + rho0 * diag(Qcell)
  }

  # ensure PSD with small jitter
  if (jitter > 0) {
    K_cell <- K_cell + jitter * diag(Qcell)
  }

  # normalization (optional)
  if (normalize != "none") {
    if (normalize == "unit_trace") {
      tr <- sum(diag(K_cell)); if (tr > 0) K_cell <- K_cell / tr
    } else if (normalize == "unit_fro") {
      nf <- sqrt(sum(K_cell^2)); if (nf > 0) K_cell <- K_cell / nf
    } else if (normalize == "max_diag") {
      md <- max(diag(K_cell)); if (md > 0) K_cell <- K_cell / md
    }
  }

  # --- 4) Map to effect-coded space if requested ---
  info <- list(levels = Ls, factor_names = fact_names, term_names = term_names,
               basis = basis, map = NULL, blocks = NULL)

  if (basis == "cell") {
    K <- K_cell
  } else {
    if (is.null(contrasts)) stop("basis='effect' requires a named list `contrasts` with L_i x d_i contrast matrices.")
    if (!all(names(contrasts) %in% fact_names)) stop("All `contrasts` names must match `factors`.")

    # stack block maps T_S for each effect block in `block_structure` (or infer from terms)
    # For main effects: T_A = C_A ⊗ 1 ⊗ ...; For A:B: T_AB = C_A ⊗ C_B ⊗ ...
    C_list <- lapply(fact_names, function(nm) {
      Ci <- contrasts[[nm]]
      if (is.null(Ci)) stop("Missing contrast for factor ", nm)
      if (!is.matrix(Ci) || nrow(Ci) != Ls[[nm]]) stop("Contrast for factor ", nm, " must be L x d.")
      Ci
    })
    names(C_list) <- fact_names

    one_vec <- function(L) matrix(1, L, 1)

    # helper to build a term map T_S of size (prod L_i) x (prod d_i for i in S)
    term_map <- function(S) {
      mats <- vector("list", length(fact_names))
      out_dim <- 1L
      for (i in seq_along(fact_names)) {
        nm <- fact_names[i]
        if (nm %in% S) {
          mats[[i]] <- C_list[[nm]]
          out_dim <- out_dim * ncol(C_list[[nm]])
        } else {
          mats[[i]] <- one_vec(Ls[[nm]])
        }
      }
      T <- kron_all(mats) # (∏ L_i) x out_dim
      attr(T, "out_dim") <- out_dim
      T
    }

    # Decide which effect blocks to include and in what order
    if (is.null(block_structure)) {
      # include any singleton S in terms and any higher-order S in terms
      block_structure <- unique(term_names)
    }
    # build and stack T blocks
    T_blocks <- lapply(terms, term_map)
    names(T_blocks) <- term_names

    # Filter/Order by block_structure
    missing <- setdiff(block_structure, names(T_blocks))
    if (length(missing) > 0) stop("Blocks requested but not present in terms: ", paste(missing, collapse=", "))
    T_list <- T_blocks[ block_structure ]

    # Build big T and block indices
    out_cols <- vapply(T_list, function(T) attr(T, "out_dim"), integer(1))
    T <- do.call(cbind, T_list)
    blocks <- split(seq_len(sum(out_cols)),
                    rep(seq_along(out_cols), times = out_cols))

    # Effect-space kernel
    K <- crossprod(T, K_cell %*% T)

    info$map <- T
    info$blocks <- blocks
  }

  list(K=K, K_cell=K_cell, info=info)
}

#' Square root and inverse square root of a PSD kernel
#' 
#' Computes the matrix square root and inverse square root of a positive 
#' semi-definite kernel using eigendecomposition with optional regularization.
#' 
#' @details
#' For a positive semi-definite matrix K, this function computes K^(1/2) and
#' K^(-1/2) using the eigendecomposition K = V * diag(lambda) * V', where
#' K^(1/2) = V * diag(sqrt(lambda)) * V' and K^(-1/2) = V * diag(1/sqrt(lambda)) * V'.
#' The jitter parameter adds a small ridge to eigenvalues to prevent numerical
#' issues with near-zero eigenvalues.
#' 
#' @param K Symmetric positive semi-definite matrix
#' @param jitter Small ridge parameter to ensure numerical stability (default 1e-10)
#' @return A list containing:
#'   \item{Khalf}{The matrix square root K^(1/2) with same dimensions as K}
#'   \item{Kihalf}{The matrix inverse square root K^(-1/2) with same dimensions as K}
#'   \item{evals}{Eigenvalues (after jitter adjustment)}
#'   \item{evecs}{Eigenvectors matrix}
#' @examples
#' # Simple 2x2 positive definite matrix
#' K <- matrix(c(2, 1, 1, 2), 2, 2)
#' result <- kernel_roots(K)
#' print(dim(result$Khalf))  # 2x2
#' @export
kernel_roots <- function(K, jitter=1e-10) {
  if (!is.matrix(K)) stop("K must be a matrix.")
  Ksym <- (K + t(K)) / 2
  ee <- eigen(Ksym, symmetric = TRUE)
  vals <- pmax(ee$values, jitter)
  V <- ee$vectors
  Khalf  <- V %*% (diag(sqrt(vals), nrow=length(vals))) %*% t(V)
  Kihalf <- V %*% (diag(1/sqrt(vals), nrow=length(vals))) %*% t(V)
  list(Khalf=Khalf, Kihalf=Kihalf, evals=vals, evecs=V)
}

#' Kernel alignment score
#' 
#' Computes the normalized Frobenius inner product between two symmetric matrices,
#' which measures the similarity between two kernels. Useful for model selection
#' and comparing different kernel parameterizations.
#' 
#' @details
#' The kernel alignment is computed as: align(A,B) = <A,B>_F / (||A||_F * ||B||_F),
#' where <A,B>_F = sum(A * B) is the Frobenius inner product and ||A||_F is the
#' Frobenius norm. This provides a normalized measure of kernel similarity.
#' 
#' @param A First symmetric matrix
#' @param B Second symmetric matrix (must have same dimensions as A)
#' @return Scalar alignment score in [-1,1], where 1 indicates perfect positive 
#'   alignment, 0 indicates orthogonality, and -1 indicates perfect negative 
#'   alignment
#' @examples
#' # Compare two 3x3 matrices
#' A <- matrix(c(2, 1, 0, 1, 2, 1, 0, 1, 2), 3, 3)
#' B <- matrix(c(1, 0, 0, 0, 1, 0, 0, 0, 1), 3, 3)
#' align_score <- kernel_alignment(A, B)
#' print(align_score)  # alignment between A and identity matrix
#' @export
kernel_alignment <- function(A, B) {
  if (!all(dim(A) == dim(B))) stop("A and B must have the same dimensions.")
  num <- sum(A * B)
  den <- sqrt(sum(A * A) * sum(B * B))
  if (den <= 0) return(0)
  num / den
}

#' Example 5x5 factorial kernel
#' 
#' Creates a simple 5x5 factorial design kernel with two nominal factors,
#' useful for testing and demonstration purposes.
#' 
#' @param rho0 Weight for identity/ridge term (default 1e-8)
#' @param rhoA Weight for main effect of factor A (default 1)
#' @param rhoB Weight for main effect of factor B (default 1)  
#' @param rhoAB Weight for A:B interaction (default 0)
#' @return A list with kernel matrices and metadata (see \code{\link{design_kernel}})
#' @examples
#' # Simple additive model (no interaction)
#' K1 <- example_kernel_5x5(rhoA=1, rhoB=1, rhoAB=0)
#' 
#' # Full factorial model
#' K2 <- example_kernel_5x5(rhoA=1, rhoB=1, rhoAB=1)
#' @export
example_kernel_5x5 <- function(rho0=1e-8, rhoA=1, rhoB=1, rhoAB=0) {
  factors <- list(A=list(L=5, type="nominal"),
                  B=list(L=5, type="nominal"))
  terms <- list("A", "B", c("A","B"))
  rho <- c("A"=rhoA, "B"=rhoB, "A:B"=rhoAB)
  design_kernel(factors, terms=terms, rho=rho, include_intercept=TRUE, rho0=rho0, basis="cell")
}

#' Build sum-to-zero contrasts
#' 
#' Creates sum-to-zero contrast matrices for a set of factors,
#' suitable for use with effect-coded kernels.
#' 
#' @param Ls Named integer vector specifying the number of levels per factor
#' @return Named list of contrast matrices, each of dimension L x (L-1) where
#'   L is the number of levels for that factor. Each matrix has row sums of zero.
#' @examples
#' # Two factors with 3 and 4 levels
#' contrasts <- sum_contrasts(c(A=3, B=4))
#' print(dim(contrasts$A))  # 3 x 2
#' print(dim(contrasts$B))  # 4 x 3
#' @export
sum_contrasts <- function(Ls) {
  lapply(Ls, function(L) {
    # contr.sum returns (L x (L-1)); ensure numeric
    cm <- contr.sum(L)
    storage.mode(cm) <- "double"
    cm
  })
}

#' Build Helmert orthonormal contrasts
#' 
#' Creates orthonormal Helmert contrast matrices for a set of factors,
#' suitable for use with effect-coded kernels. The resulting contrasts
#' are orthonormal, which can improve numerical stability.
#' 
#' @param Ls Named integer vector specifying the number of levels per factor
#' @return Named list of orthonormal contrast matrices, each of dimension 
#'   L x (L-1) where L is the number of levels for that factor. Each matrix
#'   has orthonormal columns.
#' @examples
#' # Two factors with 3 and 4 levels  
#' contrasts <- helmert_contrasts(c(A=3, B=4))
#' print(dim(contrasts$A))  # 3 x 2
#' print(dim(contrasts$B))  # 4 x 3
#' @export
helmert_contrasts <- function(Ls) {
  lapply(Ls, function(L) {
    cm <- contr.helmert(L)  # L x (L-1)
    # Orthonormalize columns (optional)
    Q <- qr.Q(qr(cm))
    storage.mode(Q) <- "double"
    Q
  })
}

