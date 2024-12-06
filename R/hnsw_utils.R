#' Perform Radius-based Search using HNSW
#'
#' This function emulates radius-based search using HNSW's k-nearest neighbors search.
#' It first finds k nearest neighbors, then filters them based on the radius threshold.
#'
#' @param query_points Matrix of query points where each row is a point
#' @param reference_points Matrix of reference points where each row is a point
#' @param radius The radius within which to find neighbors (in the same units as the distance metric)
#' @param max_neighbors Maximum number of neighbors to consider (like max_neighbour in rflann)
#' @param ef Search parameter that controls speed/accuracy trade-off (default: 100)
#' @param M Parameter that controls index size/speed trade-off (default: 16)
#' @param distance Distance metric to use: "l2" (squared Euclidean), "euclidean", "cosine", or "ip" (inner product)
#'
#' @return A list with two components:
#'   \item{indices}{A list where each element contains indices of neighbors within radius}
#'   \item{distances}{A list where each element contains distances to neighbors within radius}
#'
#' @importFrom RcppHNSW hnsw_build hnsw_search
#' @export
hnsw_radius_search <- function(query_points, reference_points, radius, 
                             max_neighbors = 100, ef = 100, M = 16,
                             distance = "l2") {
    
    # Build the index
    ann <- hnsw_build(reference_points, distance = distance, M = M, ef = ef)
    
    # For L2 distance, we need to square the radius as HNSW uses squared Euclidean
    search_radius <- if(distance == "l2") radius^2 else radius
    
    # Search for k nearest neighbors
    nn_results <- hnsw_search(query_points, ann, k = max_neighbors)
    
    # Filter results by radius
    indices <- lapply(1:nrow(query_points), function(i) {
        idx <- nn_results$idx[i,]
        dist <- nn_results$dist[i,]
        keep <- dist <= search_radius
        idx[keep]
    })
    
    distances <- lapply(1:nrow(query_points), function(i) {
        dist <- nn_results$dist[i,]
        keep <- dist <= search_radius
        dist[keep]
    })
    
    list(indices = indices, distances = distances)
}
