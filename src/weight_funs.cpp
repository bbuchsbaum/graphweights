// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
#include <unordered_set>
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

// [[Rcpp::export]]
List expand_similarity_cpp(IntegerVector indices, NumericMatrix simmat, double thresh ) {
  int n = indices.length();
  //List triplets;

  std::vector<NumericVector> triplets;

  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      if (j > i) {
        continue;
      } else {
        int r = indices(i)-1;
        int c = indices(j)-1;
        if (r < 0 || c < 0) {
          continue;
        }
        if (simmat(r,c) > thresh) {
          NumericVector ii = NumericVector::create(i+1,j+1,simmat(r,c));
          triplets.push_back(ii);
        }

      }
    }
  }
  return Rcpp::wrap(triplets);
}

// [[Rcpp::export]]
List expand_similarity_below_cpp(IntegerVector indices, NumericMatrix simmat, double thresh ) {
  int n = indices.length();
  //List triplets;

  std::vector<NumericVector> triplets;

  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      if (j > i) {
        continue;
      } else {
        int r = indices(i)-1;
        int c = indices(j)-1;
        if (simmat(r,c) < thresh) {
          NumericVector ii = NumericVector::create(i+1,j+1,simmat(r,c));
          triplets.push_back(ii);
        }

      }
    }
  }
  return Rcpp::wrap(triplets);
}


// [[Rcpp::export]]
double norm_heat_kernel(NumericVector x1, NumericVector x2, double sigma) {
  double norm_dist = sqrt(sum(pow((x1-x2),2))/(2*x1.length()));
  return exp(-norm_dist/(2*pow(sigma,2)));
}



NumericVector distance_heat(NumericVector dist, double sigma) {
  return exp(-dist/(2*pow(sigma,2)));
}



NumericVector normalized_distance(NumericVector dist, double sigma, int len) {
  NumericVector norm_dist = pow(dist,2)/(2*len);
  return exp(-norm_dist/(2*pow(sigma,2)));
}

// [[Rcpp::export]]
IntegerVector order_vec(NumericVector x) {
  NumericVector sorted = clone(x).sort();
  return match(sorted, x);
}

// [[Rcpp::export]]
NumericMatrix cross_fspatial_weights(List indices, List distances, NumericMatrix feature_mat1,
                                     NumericMatrix feature_mat2,
                                     double nels, double sigma,
                                     double fsigma, double alpha,
                                     int maxk,
                                     bool binary) {

  int n = indices.length();
  List out(n);

  NumericMatrix wout(nels, 3);

  int count = 0;

  double alpha2 = 1 - alpha;

  for(int i = 0; i < n; ++i) {
    SEXP ll = indices[i];
    IntegerVector ind(ll);

    if (ind.size() > 0) {
      NumericVector vals;
      NumericVector fvals(ind.size());

      if (binary) {
        vals = rep(1.0, ind.size());
      } else {
        vals = distance_heat(distances[i], sigma);
      }

      NumericMatrix::Row f1 = feature_mat1(i,_);
      if (maxk < ind.size()) {
        for (int j=0; j<vals.length(); j++) {
          NumericMatrix::Row f2 = feature_mat2(ind[j]-1,_);
          double s = norm_heat_kernel(f1, f2, fsigma);
          fvals[j] = alpha*vals[j] + alpha2*s;
        }

        IntegerVector ord = order_vec(fvals);
        for (int j=0; j<maxk; j++) {
          wout(count, 0) = i+1;
          wout(count, 1) = ind[ord[j]-1];
          wout(count, 2) = fvals[ord[j]-1];
          count++;
        }
      } else {

        for (int j=0; j<vals.length(); j++) {
          NumericMatrix::Row f2 = feature_mat2(ind[j]-1,_);
          double s = norm_heat_kernel(f1, f2, fsigma);
          wout(count, 0) = i+1;
          wout(count, 1) = ind[j];
          wout(count, 2) = alpha*vals[j] + alpha2*s;
          count++;
        }
      }
    }
  }

  return wout;

}


// [[Rcpp::export]]
NumericMatrix bilateral_weights(List indices, List distances, NumericMatrix feature_mat,
                               double nels, double sigma, double fsigma) {
  int n = indices.length();
  List out(n);

  NumericMatrix wout(nels, 3);

  int count = 0;


  for(int i = 0; i < n; ++i) {
    SEXP ll = indices[i];
    IntegerVector ind(ll);
    // first compute spatial similarity
    if (ind.size() > 0) {
      NumericVector vals;

      vals = distance_heat(distances[i], sigma);

      NumericMatrix::Row f1 = feature_mat(i,_);

      for (int j=0; j<vals.length(); j++) {
        NumericMatrix::Row f2 = feature_mat(ind[j]-1,_);
        double s = norm_heat_kernel(f1, f2, fsigma);
        wout(count, 0) = i+1;
        wout(count, 1) = ind[j];
        wout(count, 2) = vals[j]*s;
        count++;
      }
    }
  }

  return wout;

}


// [[Rcpp::export]]
NumericMatrix fspatial_weights(List indices, List distances, NumericMatrix feature_mat,
                               double nels, double sigma, double fsigma, double alpha, bool binary) {
  int n = indices.length();
  List out(n);

  NumericMatrix wout(nels, 3);

  int count = 0;

  double alpha2 = 1 - alpha;

  for(int i = 0; i < n; ++i) {
    SEXP ll = indices[i];
    IntegerVector ind(ll);
    // first compute spatial similarity
    if (ind.size() > 0) {
      NumericVector vals;
      if (binary) {
        vals = rep(1.0, ind.size());
      } else {
        vals = distance_heat(distances[i], sigma);
      }

      NumericMatrix::Row f1 = feature_mat(i,_);

      for (int j=0; j<vals.length(); j++) {
        NumericMatrix::Row f2 = feature_mat(ind[j]-1,_);
        double s = norm_heat_kernel(f1, f2, fsigma);
        wout(count, 0) = i+1;
        wout(count, 1) = ind[j];
        wout(count, 2) = alpha*vals[j] + alpha2*s;
        count++;
      }
    }
  }

  return wout;

}

// [[Rcpp::export]]
NumericMatrix spatial_weights(List indices, List distances, double nels, double sigma, bool binary) {
  int n = indices.length();
  List out(n);

  NumericMatrix wout(nels, 3);

  int count = 0;

  for(int i = 0; i < n; ++i) {
    SEXP ll = indices[i];
    IntegerVector ind(ll);

    if (ind.size() > 0) {
      // compute coordinate distances
      NumericVector vals;
      if (binary) {
        vals = rep(1.0, ind.size());
      } else {
        vals = distance_heat(distances[i], sigma);
      }

      for (int j=0; j<ind.length(); j++) {
          wout(count, 0) = i+1;
          wout(count, 1) = ind[j];
          wout(count, 2) = vals[j];
          count++;
      }
    }
  }

  return wout;
}

