#' @import RcppHNSW
nnsearcher <- function(X, labels=1:nrow(labels), ...,
                       distance=c("l2", "euclidean", "cosine", "ip"), M=16, ef=200) {
  distance <- match.arg(distance)
  X <- as.matrix(X)
  ann <- RcppHNSW::hnsw_build(X, distance, ef=ef, M=M)

  structure(list(
    X=X,
    labels=labels,
    ann=ann,
    distance=distance,
    ef=ef,
    M=M),
    class="nnsearcher")
}

#' @export
search_result.nnsearcher <- function(x, result) {
  attr(result, "len") <- ncol(x$X)
  attr(result, "metric") <- x$distance
  class(result) <- "nn_search"
  result
}

#' @export
dist_to_sim.nn_search <- function(x, method = c("heat", "binary", "normalized", "cosine", "correlation"), sigma=1) {
  method <- match.arg(method)
  fun <- get_neighbor_fun(method, x$len, x$sigma)
  x$dist <- fun(x$dist)
  x
}


#' @export
dist_to_sim.Matrix <- function(x, method = c("heat", "binary", "normalized", "cosine", "correlation"), sigma=1, len=1) {
  method <- match.arg(method)
  fun <- get_neighbor_fun(method, len, sigma)

  wh <- which(x != 0)
  v <- fun(x[wh])
  x[wh] <- v
  x
}

#' @export
adjacency.nnsearch <- function(x, idim=nrow(x$idx), jdim=max(x$idx), return_triplet=FALSE) {
  indices_to_sparse(as.matrix(x$idx), as.matrix(x$dist), return_triplet=FALSE, idim=idim, jdim=jdim)
}



#' @export
find_nn.nnsearcher <- function(x, query=NULL, k=5) {
  ret <- if (!is.null(query)) {
    chk::chk_matrix(query)
    RcppHNSW::hnsw_search(query, x$ann, k = k)
  } else {
    RcppHNSW::hnsw_search(x$X, x$ann, k = k)
  }

  ret$labels <- t(apply(ret$idx, 1,  function(i) x$labels[i]))
  search_result(x,ret)
}

#' @export
find_nn_among.nn_searcher <- function(x, k=5, idx) {
  chk::chk_numeric(idx)

  X1 <- x$X[idx,]
  ann <- RcppHNSW::hnsw_build(X1, x$distance, M=x$M, ef=x$ef)
  nnres <- RcppHNSW::hnsw_search(X1, ann, k=k)
  class(nnres) <- "nn_search"
  nnres
}

find_nn_among.class_graph <- function(x, X, k=5) {
  searcher <- nnsearcher(X,x$labels)
  ret <- lapply(x$class_indices, function(ind) {
    nnr <- find_nn_among.nn_searcher(searcher, k=k, ind)
    nnr$idx <- do.call(rbind, lapply(1:nrow(nnr$idx), function(i) {
      ind[nnr$idx[i,]]
    }))

    nnr$labels <- t(apply(nnr$idx, 1,  function(i) x$labels[i]))
    nnr
  })

  idx <- do.call(rbind, lapply(ret, "[[", "idx"))
  d <- do.call(rbind, lapply(ret, "[[", "dist"))
  labels <- do.call(rbind, lapply(ret, "[[", "labels"))
  search_result(searcher, list(idx=idx, dist=d, labels=labels))
}

#' @export
find_nn_between.nn_searcher <- function(x, k=5, idx1, idx2, restricted=FALSE) {
  chk::chk_numeric(idx1)
  chk::chk_numeric(idx2)

  ret <- if (!restricted) {
    X1 <- x$X[idx1,]
    X2 <- x$X[idx2,]
    ann <- RcppHNSW::hnsw_build(X1, x$distance, M=x$M, ef=x$ef)
    nnres <- RcppHNSW::hnsw_search(X2, ann, k=k)
    nnres$idx <- do.call(rbind, lapply(1:nrow(nnres$idx), function(i) {
      idx1[nnres$idx[i,]]
    }))
    nnres$labels <- t(apply(nnres$idx, 1,  function(i) x$labels[i]))
    nnres
  } else {
    find_nn.nnsearcher(x, x$X[idx2,])
  }

  class(nnres) <- "nn_search"
  nnres
}

#' @export
neighbor_graph.nn_searcher <- function(x, query=NULL, k=5, type=c("normal", "asym", "mutual"),
                                                           transform=c("heat", "binary", "euclidean",
                                                              "normalized", "cosine", "correlation"), sigma=1) {

  nn <- find_nn.nnsearcher(x,k=(k+1))
  D <- nn$dist[,2:(k+1)]

  nnd <- sqrt(nn$dist[, 2:ncol(nn$dist),drop=FALSE])
  nni <- nn$idx[, 2:ncol(nn$idx), drop=FALSE]


  hfun <- get_neighbor_fun(transform, len=ncol(x$X), sigma=sigma)
  hval <- hfun(D)

  W <- if (k == 1) {
    indices_to_sparse(as.matrix(nni), as.matrix(hval), idim=nrow(X), jdim=nrow(X))
  } else {
    indices_to_sparse(nni, hval, idim=nrow(X), jdim=nrow(X))
  }

  gg <- if (type == "normal") {
    igraph::graph_from_adjacency_matrix(W, weighted=TRUE, mode="max")
  } else if (type == "mutual") {
    igraph::graph_from_adjacency_matrix(W, weighted=TRUE, mode="min")
  } else if (type == "asym") {
    igraph::graph_from_adjacency_matrix(W, weighted=TRUE, mode="directed")
  }

  neighbor_graph(gg, params=list(k=k,
                                transform=transform,
                                sigma=sigma,
                                type=type,
                                labels=x$labels))


}
