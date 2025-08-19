// [[Rcpp::plugins(cpp11)]]
#include <RcppArmadillo.h>
#include <vector>
#include <utility> // For std::pair
#include <algorithm> // For std::sort
#include <cmath> // For std::sqrt, std::exp
#include <limits> // For std::numeric_limits
#include <numeric> // For std::iota

using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

//
//double heat_kernel(NumericVector x1, NumericVector x2, double sigma) {
//  double dist = sqrt(sum(pow((x1-x2),2)));
//  return exp(-dist / (2*(pow(sigma,2))));
//}

// Optimized distance_heat
inline NumericVector distance_heat(NumericVector dist, double sigma) {
  int n = dist.size();
  NumericVector out(n);
  double denom = 2.0 * sigma * sigma;
  if (denom == 0) { // Avoid division by zero
      stop("sigma cannot be zero in distance_heat");
  }
  for (int i = 0; i < n; ++i) {
    out[i] = std::exp(-dist[i] / denom);
  }
  return out;
}

// Optimized norm_heat_kernel
inline double norm_heat_kernel(const NumericVector& x1, const NumericVector& x2, double fsigma) {
  int len = x1.size();
  if (len != x2.size() || len == 0) {
      stop("Vectors must have the same non-zero length in norm_heat_kernel");
  }
  double dist_sq = 0.0;
  for (int i = 0; i < len; ++i) {
    double diff = x1[i] - x2[i];
    dist_sq += diff * diff;
  }
  // Note: Original formula had dist = sqrt(sum(pow(...))/(2*len)).
  // Here we implement sqrt(dist_sq / (2.0 * len))
  double norm_dist = std::sqrt(dist_sq / (2.0 * len));
  double denom = 2.0 * fsigma * fsigma;
   if (denom == 0) { // Avoid division by zero
      stop("fsigma cannot be zero in norm_heat_kernel");
  }
  return std::exp(-norm_dist / denom);
}

// [[Rcpp::export]]
NumericMatrix expand_similarity_cpp(IntegerVector indices, NumericMatrix simmat, double thresh) {
  int n = indices.size();
  int N_sim = simmat.nrow();
  std::vector<std::vector<double>> triplets;

  for (int i = 0; i < n; ++i) {
    int r = indices[i] - 1; // Convert to 0-based index
    if (r < 0 || r >= N_sim) continue; // Bounds check

    for (int j = 0; j <= i; ++j) { // Iterate only lower triangle including diagonal
      int c = indices[j] - 1; // Convert to 0-based index
      if (c < 0 || c >= N_sim) continue; // Bounds check

      // Ensure we access within matrix bounds if simmat isn't square or fully populated
      if (r >= simmat.nrow() || c >= simmat.ncol()) continue;

      if (simmat(r, c) > thresh) {
        triplets.push_back({(double)(i + 1), (double)(j + 1), simmat(r, c)}); // Store 1-based indices
      }
    }
  }

  int m = triplets.size();
  NumericMatrix out(m, 3);
  for (int i = 0; i < m; ++i) {
    out(i, 0) = triplets[i][0];
    out(i, 1) = triplets[i][1];
    out(i, 2) = triplets[i][2];
  }

  return out;
}

// [[Rcpp::export]]
NumericMatrix expand_similarity_below_cpp(IntegerVector indices, NumericMatrix simmat, double thresh) {
  int n = indices.size();
  int N_sim = simmat.nrow();
  std::vector<std::vector<double>> triplets;

  for (int i = 0; i < n; ++i) {
    int r = indices[i] - 1; // Convert to 0-based index
    if (r < 0 || r >= N_sim) continue; // Bounds check

    for (int j = 0; j <= i; ++j) { // Iterate only lower triangle including diagonal
      int c = indices[j] - 1; // Convert to 0-based index
      if (c < 0 || c >= N_sim) continue; // Bounds check

       // Ensure we access within matrix bounds if simmat isn't square or fully populated
      if (r >= simmat.nrow() || c >= simmat.ncol()) continue;

      if (simmat(r, c) < thresh) { // Changed condition to '<'
        triplets.push_back({(double)(i + 1), (double)(j + 1), simmat(r, c)}); // Store 1-based indices
      }
    }
  }

  int m = triplets.size();
  NumericMatrix out(m, 3);
  for (int i = 0; i < m; ++i) {
    out(i, 0) = triplets[i][0];
    out(i, 1) = triplets[i][1];
    out(i, 2) = triplets[i][2];
  }

  return out;
}


// Fixed order_vec using std::sort
// [[Rcpp::export]]
IntegerVector order_vec(NumericVector x) {
  int n = x.size();
  // Create vector of pairs (value, original_index)
  std::vector<std::pair<double, int>> val_idx(n);
  for (int i = 0; i < n; ++i) {
    val_idx[i] = std::make_pair(x[i], i);
  }

  // Sort based on values (first element of pair)
  std::sort(val_idx.begin(), val_idx.end());

  // Extract the sorted original indices
  IntegerVector out(n);
  for (int i = 0; i < n; ++i) {
    out[i] = val_idx[i].second + 1; // Return 1-based indices
  }
  return out;
}

// Structure to hold triplet data
struct WeightTriplet {
    double i, j, weight;
};

// [[Rcpp::export]]
NumericMatrix cross_fspatial_weights(List indices, List distances, NumericMatrix feature_mat1,
                                     NumericMatrix feature_mat2,
                                     double sigma,
                                     double fsigma, double alpha,
                                     int maxk,
                                     bool binary) {

  int n = indices.size();
  if (n != feature_mat1.nrow()) {
      stop("Length of indices/distances must match number of rows in feature_mat1");
  }

  // Use std::vector to store results temporarily
  std::vector<WeightTriplet> triplets;
  // Reserve space based on expectation (nels was used before)
  // A rough guess might be n * avg_neighbors (e.g., n * maxk if maxk is small)
  // Or start smaller and let it grow.
  triplets.reserve(n * (maxk > 0 ? maxk : 10)); // Example reservation

  double alpha2 = 1.0 - alpha;
  int N_feat2 = feature_mat2.nrow();

  for (int i = 0; i < n; ++i) {
    IntegerVector ind = Rcpp::as<IntegerVector>(indices[i]);
    NumericVector dist = Rcpp::as<NumericVector>(distances[i]);
    int k_neighbors = ind.size();

    if (k_neighbors == 0 || k_neighbors != dist.size()) continue; // Skip if no neighbors or mismatched sizes

    NumericVector spatial_vals = binary ? NumericVector(k_neighbors, 1.0) : distance_heat(dist, sigma);
    NumericVector f1 = feature_mat1(i, _);

    std::vector<std::pair<double, int>> neighbor_data; // Store {combined_value, original_index_in_ind}
    if (maxk > 0 && maxk < k_neighbors) {
        neighbor_data.reserve(k_neighbors);
    }

    for (int j = 0; j < k_neighbors; ++j) {
      int neighbor_idx = ind[j] - 1; // 0-based index for feature_mat2

      // Bounds check for feature_mat2
      if (neighbor_idx < 0 || neighbor_idx >= N_feat2) continue;

      NumericVector f2 = feature_mat2(neighbor_idx, _);
      double feature_sim = norm_heat_kernel(f1, f2, fsigma);
      double combined_val = alpha * spatial_vals[j] + alpha2 * feature_sim;

      if (R_finite(combined_val)) { // Check for NaN/Inf
          if (maxk > 0 && maxk < k_neighbors) {
              neighbor_data.push_back({combined_val, j});
          } else { // maxk is non-limiting or disabled
              triplets.push_back({(double)(i + 1), (double)ind[j], combined_val});
          }
      } else {
           warning("Non-finite weight encountered for pair (%d, %d)", i + 1, ind[j]);
      }
    }

    // If maxk is active, sort and select top k
    if (maxk > 0 && maxk < k_neighbors && !neighbor_data.empty()) {
        // Sort by combined_val (descending)
        std::sort(neighbor_data.rbegin(), neighbor_data.rend()); // Sort descending by weight

        int num_to_take = std::min((int)neighbor_data.size(), maxk);
        for (int j = 0; j < num_to_take; ++j) {
            int original_j = neighbor_data[j].second;
            triplets.push_back({(double)(i + 1), (double)ind[original_j], neighbor_data[j].first});
        }
    }
  }

  // Create the final matrix from the collected triplets
  int m = triplets.size();
  NumericMatrix wout(m, 3);
  for(int i = 0; i < m; ++i) {
      wout(i, 0) = triplets[i].i;
      wout(i, 1) = triplets[i].j;
      wout(i, 2) = triplets[i].weight;
  }
  return wout;
}


// [[Rcpp::export]]
NumericMatrix bilateral_weights(List indices, List distances, NumericMatrix feature_mat,
                                double sigma, double fsigma) {
  int n = indices.size();
  if (n != feature_mat.nrow()) {
      stop("Length of indices/distances must match number of rows in feature_mat");
  }
  int N_feat = feature_mat.nrow();

  std::vector<WeightTriplet> triplets;
  // Guess initial size - can adjust
  triplets.reserve(n * 10); 

  for (int i = 0; i < n; ++i) {
    IntegerVector ind = Rcpp::as<IntegerVector>(indices[i]);
    NumericVector dist = Rcpp::as<NumericVector>(distances[i]);
    int k_neighbors = ind.size();

    if (k_neighbors == 0 || k_neighbors != dist.size()) continue;

    NumericVector spatial_vals = distance_heat(dist, sigma);
    NumericVector f1 = feature_mat(i, _);

    for (int j = 0; j < k_neighbors; ++j) {
      int neighbor_idx = ind[j] - 1; // 0-based index

      // Bounds check
      if (neighbor_idx < 0 || neighbor_idx >= N_feat) continue;

      NumericVector f2 = feature_mat(neighbor_idx, _);
      double feature_sim = norm_heat_kernel(f1, f2, fsigma);
      double final_weight = spatial_vals[j] * feature_sim;

      if (R_finite(final_weight)) {
          triplets.push_back({(double)(i + 1), (double)ind[j], final_weight});
      } else {
           warning("Non-finite weight encountered for pair (%d, %d)", i + 1, ind[j]);
      }
    }
  }

  // Create the final matrix
  int m = triplets.size();
  NumericMatrix wout(m, 3);
  for(int i = 0; i < m; ++i) {
      wout(i, 0) = triplets[i].i;
      wout(i, 1) = triplets[i].j;
      wout(i, 2) = triplets[i].weight;
  }
  return wout;
}


// [[Rcpp::export]]
NumericMatrix fspatial_weights(List indices, List distances, NumericMatrix feature_mat,
                               double sigma, double fsigma, double alpha, bool binary) {
  int n = indices.size();
   if (n != feature_mat.nrow()) {
      stop("Length of indices/distances must match number of rows in feature_mat");
  }
  int N_feat = feature_mat.nrow();

  std::vector<WeightTriplet> triplets;
  triplets.reserve(n * 10); 

  double alpha2 = 1.0 - alpha;

  for (int i = 0; i < n; ++i) {
    IntegerVector ind = Rcpp::as<IntegerVector>(indices[i]);
    NumericVector dist = Rcpp::as<NumericVector>(distances[i]);
    int k_neighbors = ind.size();

    if (k_neighbors == 0 || k_neighbors != dist.size()) continue;

    NumericVector spatial_vals = binary ? NumericVector(k_neighbors, 1.0) : distance_heat(dist, sigma);
    NumericVector f1 = feature_mat(i, _);

    for (int j = 0; j < k_neighbors; ++j) {
      int neighbor_idx = ind[j] - 1; // 0-based index

      // Bounds check
      if (neighbor_idx < 0 || neighbor_idx >= N_feat) continue;

      NumericVector f2 = feature_mat(neighbor_idx, _);
      double feature_sim = norm_heat_kernel(f1, f2, fsigma);
      double final_weight = alpha * spatial_vals[j] + alpha2 * feature_sim;

       if (R_finite(final_weight)) {
          triplets.push_back({(double)(i + 1), (double)ind[j], final_weight});
      } else {
          warning("Non-finite weight encountered for pair (%d, %d)", i + 1, ind[j]);
      }
    }
  }

  // Create the final matrix
  int m = triplets.size();
  NumericMatrix wout(m, 3);
  for(int i = 0; i < m; ++i) {
      wout(i, 0) = triplets[i].i;
      wout(i, 1) = triplets[i].j;
      wout(i, 2) = triplets[i].weight;
  }
  return wout;
}

// [[Rcpp::export]]
NumericMatrix spatial_weights(List indices, List distances, double sigma, bool binary) {
  int n = indices.size();
  
  std::vector<WeightTriplet> triplets;
  triplets.reserve(n * 10);

  for (int i = 0; i < n; ++i) {
    IntegerVector ind = Rcpp::as<IntegerVector>(indices[i]);
    NumericVector dist = Rcpp::as<NumericVector>(distances[i]);
    int k_neighbors = ind.size();

    if (k_neighbors == 0 || k_neighbors != dist.size()) continue;

    NumericVector spatial_vals = binary ? NumericVector(k_neighbors, 1.0) : distance_heat(dist, sigma);

    for (int j = 0; j < k_neighbors; ++j) {
      // No feature matrix index check needed here
      double final_weight = spatial_vals[j];

       if (R_finite(final_weight)) {
            triplets.push_back({(double)(i + 1), (double)ind[j], final_weight});
       } else {
            warning("Non-finite weight encountered for pair (%d, %d)", i + 1, ind[j]);
       }
    }
  }

  // Create the final matrix
  int m = triplets.size();
  NumericMatrix wout(m, 3);
  for(int i = 0; i < m; ++i) {
      wout(i, 0) = triplets[i].i;
      wout(i, 1) = triplets[i].j;
      wout(i, 2) = triplets[i].weight;
  }
  return wout;
}

