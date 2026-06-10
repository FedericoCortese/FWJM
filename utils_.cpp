#include <Rcpp.h>
#include <algorithm>
#include <vector>
#include <cmath>
#include <numeric>

using namespace Rcpp;

// helper to compute median of a std::vector<double>
double median_vec(std::vector<double>& v) {
  std::sort(v.begin(), v.end());
  int n = v.size();
  if (n % 2 == 1) {
    return v[(n - 1) / 2];
  } else {
    return 0.5 * (v[n/2 - 1] + v[n/2]);
  }
}

void median_halves_q1q3(const std::vector<double>& vals, double &q1, double &q3);

Rcpp::NumericMatrix gower_dist(const Rcpp::NumericMatrix& Y,
                               const Rcpp::NumericMatrix& mu,
                               Rcpp::Nullable<Rcpp::IntegerVector> feat_type = R_NilValue,
                               std::string scale = "m");

// restituisce Q1 e Q3 usando median-halves
void median_halves_q1q3(const std::vector<double>& vals, double &q1, double &q3) {
  std::vector<double> v = vals;
  std::sort(v.begin(), v.end());
  int n = v.size();
  
  if (n == 0) {
    q1 = R_NaReal;
    q3 = R_NaReal;
    return;
  }
  
  int half = n / 2;
  std::vector<double> lower, upper;
  
  if (n % 2 == 0) {
    lower.assign(v.begin(), v.begin() + half);
    upper.assign(v.begin() + half, v.end());
  } else {
    lower.assign(v.begin(), v.begin() + half);
    upper.assign(v.begin() + half + 1, v.end());
  }
  
  if (lower.empty()) q1 = R_NaReal; else q1 = median_vec(lower);
  if (upper.empty()) q3 = R_NaReal; else q3 = median_vec(upper);
}

// [[Rcpp::export]]
NumericMatrix gower_dist(const NumericMatrix& Y,
                         const NumericMatrix& mu,
                         Nullable<IntegerVector> feat_type,
                         std::string scale ) {
  int n = Y.nrow();
  int p = Y.ncol();
  int m = mu.nrow();
  
  IntegerVector ft;
  if (feat_type.isNotNull()) {
    ft = feat_type.get();
    if (ft.size() != p)
      stop("feat_type must have length = ncol(Y)");
  } else {
    ft = IntegerVector(p, 0);
  }
  
  std::vector<double> s_p(p);
  for (int j = 0; j < p; ++j) {
    if (ft[j] == 0) {
      std::vector<double> col(n);
      for (int i = 0; i < n; ++i) col[i] = Y(i, j);
      
      if (scale == "m") {
        double mn = col[0], mx = col[0];
        for (double v : col) {
          if (v < mn) mn = v;
          if (v > mx) mx = v;
        }
        s_p[j] = mx - mn;
      } else if (scale == "i") {
        std::sort(col.begin(), col.end());
        double q1 = col[(int)std::floor((n - 1) * 0.25)];
        double q3 = col[(int)std::floor((n - 1) * 0.75)];
        s_p[j] = (q3 - q1) / 1.35;
      } else if (scale == "s") {
        double sum = 0.0;
        for (double v : col) sum += v;
        double mu_col = sum / n;
        double ss = 0.0;
        for (double v : col) ss += (v - mu_col) * (v - mu_col);
        s_p[j] = std::sqrt(ss / (n - 1));
      } else {
        stop("Invalid scale flag: must be 'm', 'i', or 's'");
      }
      
      if (s_p[j] == 0.0) s_p[j] = 1.0;
    } else {
      s_p[j] = 1.0;
    }
  }
  
  std::vector< std::vector<int> > ord_rank_Y(p);
  std::vector< std::vector<int> > ord_rank_mu(p);
  std::vector<int> M(p, 1);
  
  for (int j = 0; j < p; ++j) {
    if (ft[j] == 2) {
      std::vector<double> vals;
      vals.reserve(n + m);
      for (int i = 0; i < n; ++i) vals.push_back(Y(i, j));
      for (int u = 0; u < m; ++u) vals.push_back(mu(u, j));
      
      std::sort(vals.begin(), vals.end());
      vals.erase(std::unique(vals.begin(), vals.end()), vals.end());
      
      int levels = vals.size();
      M[j] = levels > 1 ? levels : 1;
      
      ord_rank_Y[j].resize(n);
      ord_rank_mu[j].resize(m);
      
      for (int i = 0; i < n; ++i) {
        ord_rank_Y[j][i] =
          std::lower_bound(vals.begin(), vals.end(), Y(i, j)) - vals.begin();
      }
      for (int u = 0; u < m; ++u) {
        ord_rank_mu[j][u] =
          std::lower_bound(vals.begin(), vals.end(), mu(u, j)) - vals.begin();
      }
    }
  }
  
  NumericMatrix V(n, m);
  for (int i = 0; i < n; ++i) {
    for (int u = 0; u < m; ++u) {
      double acc = 0.0;
      for (int j = 0; j < p; ++j) {
        double diff;
        if (ft[j] == 0) {
          diff = std::abs(Y(i, j) - mu(u, j)) / s_p[j];
        } else if (ft[j] == 1) {
          diff = (Y(i, j) != mu(u, j)) ? 1.0 : 0.0;
        } else {
          double denom = double(M[j] - 1);
          diff = denom > 0.0 ? std::abs(ord_rank_Y[j][i] - ord_rank_mu[j][u]) / denom : 0.0;
        }
        acc += diff;
      }
      V(i, u) = acc / p;
    }
  }
  
  return V;
}

// [[Rcpp::export]]
IntegerVector initialize_states(const NumericMatrix& Y,
                                int K,
                                Nullable<IntegerVector> feat_type,
                                int reps,
                                std::string scale ) {
  int TT = Y.nrow();
  int P  = Y.ncol();
  
  IntegerVector ft;
  if (feat_type.isNotNull()) {
    ft = feat_type.get();
    if ((int)ft.size() != P) stop("feat_type must have length = ncol(Y)");
  } else {
    ft = IntegerVector(P, 0);
  }
  
  NumericMatrix Dall = gower_dist(Y, Y, ft, scale);
  
  double best_sum = R_PosInf;
  IntegerVector best_assign(TT);
  
  for (int rep = 0; rep < reps; ++rep) {
    std::vector<int> centIdx;
    centIdx.reserve(K);
    
    int idx0 = std::floor(R::runif(0, TT));
    centIdx.push_back(idx0);
    
    std::vector<double> closestDist(TT);
    for (int j = 0; j < TT; ++j) closestDist[j] = Dall(idx0, j);
    
    for (int k = 1; k < K; ++k) {
      double sumd = std::accumulate(closestDist.begin(), closestDist.end(), 0.0);
      
      if (sumd <= 0) {
        idx0 = std::floor(R::runif(0, TT));
      } else {
        double u = R::runif(0, sumd);
        double cum = 0;
        int idx = 0;
        for (; idx < TT; ++idx) {
          cum += closestDist[idx];
          if (cum >= u) break;
        }
        if (idx >= TT) idx = TT - 1;
        idx0 = idx;
      }
      
      centIdx.push_back(idx0);
      
      for (int j = 0; j < TT; ++j) {
        closestDist[j] = std::min(closestDist[j], Dall(centIdx[k], j));
      }
    }
    
    double sum_intra = 0.0;
    IntegerVector assign(TT);
    
    for (int i = 0; i < TT; ++i) {
      int best_k = 0;
      double best_d = Dall(centIdx[0], i);
      
      for (int k = 1; k < K; ++k) {
        double d = Dall(centIdx[k], i);
        if (d < best_d) {
          best_d = d;
          best_k = k;
        }
      }
      
      assign[i] = best_k + 1;
      sum_intra += best_d;
    }
    
    if (sum_intra < best_sum) {
      best_sum = sum_intra;
      best_assign = assign;
    }
  }
  
  return best_assign;
}

// [[Rcpp::export]]
Rcpp::NumericVector tukey_biw_vec_cpp(Rcpp::NumericVector u, double c = 4.685) {
  int n = u.size();
  Rcpp::NumericVector out(n);
  double c2_over6 = (c * c) / 6.0;
  
  for (int i = 0; i < n; ++i) {
    if (Rcpp::NumericVector::is_na(u[i])) {
      out[i] = NA_REAL;
    } else {
      double v = std::abs(u[i]);
      if (v <= c) {
        double w = v / c;
        double tmp = 1.0 - w * w;
        out[i] = c2_over6 * (1.0 - tmp * tmp * tmp);
      } else {
        out[i] = c2_over6;
      }
    }
  }
  
  return out;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix safe_scale_slice_median_cpp(Rcpp::NumericMatrix mat, std::string scale = "i") {
  int nr = mat.nrow();
  int nc = mat.ncol();
  Rcpp::NumericMatrix out = Rcpp::clone(mat);
  
  if (scale == "m") {
    double minv = R_PosInf;
    double maxv = R_NegInf;
    bool any_non_na = false;
    
    for (int i = 0; i < nr; ++i) {
      for (int j = 0; j < nc; ++j) {
        double val = mat(i, j);
        if (!Rcpp::NumericVector::is_na(val)) {
          any_non_na = true;
          if (val < minv) minv = val;
          if (val > maxv) maxv = val;
        }
      }
    }
    
    if (!any_non_na) return out;
    
    double sc = maxv - minv;
    if (sc > 0 && R_finite(sc)) {
      for (int i = 0; i < nr; ++i) {
        for (int j = 0; j < nc; ++j) {
          double val = out(i, j);
          if (!Rcpp::NumericVector::is_na(val)) out(i, j) = val / sc;
        }
      }
    }
    
    return out;
  }
  
  if (scale == "i") {
    std::vector<double> vals;
    vals.reserve(nr * nc);
    
    for (int i = 0; i < nr; ++i) {
      for (int j = 0; j < nc; ++j) {
        double val = mat(i, j);
        if (!Rcpp::NumericVector::is_na(val)) vals.push_back(val);
      }
    }
    
    if (vals.empty()) return out;
    
    double q1 = NA_REAL, q3 = NA_REAL;
    median_halves_q1q3(vals, q1, q3);
    if (Rcpp::NumericVector::is_na(q1) || Rcpp::NumericVector::is_na(q3)) return out;
    
    double iq = q3 - q1;
    if (iq > 0 && R_finite(iq)) {
      double scale_factor = 1.35 / iq;
      for (int i = 0; i < nr; ++i) {
        for (int j = 0; j < nc; ++j) {
          double val = out(i, j);
          if (!Rcpp::NumericVector::is_na(val)) out(i, j) = val * scale_factor;
        }
      }
    }
    
    return out;
  }
  
  return out;
}

// [[Rcpp::export]]
NumericVector gower_dist_array(const NumericMatrix& Y,
                               const NumericMatrix& mu,
                               Nullable<IntegerVector> feat_type = R_NilValue,
                               std::string scale = "m") {
  int n = Y.nrow();
  int p = Y.ncol();
  int m = mu.nrow();
  
  IntegerVector ft;
  if (feat_type.isNotNull()) {
    ft = feat_type.get();
    if (ft.size() != p) stop("feat_type must have length = ncol(Y)");
  } else {
    ft = IntegerVector(p, 0);
  }
  
  std::vector<double> s_p(p);
  for (int j = 0; j < p; ++j) {
    if (ft[j] == 0) {
      std::vector<double> col(n);
      for (int i = 0; i < n; ++i) col[i] = Y(i, j);
      
      if (scale == "m") {
        double mn = col[0], mx = col[0];
        for (double v : col) {
          if (v < mn) mn = v;
          if (v > mx) mx = v;
        }
        s_p[j] = mx - mn;
      } else if (scale == "i") {
        std::sort(col.begin(), col.end());
        double q1 = col[(int)std::floor((n - 1) * 0.25)];
        double q3 = col[(int)std::floor((n - 1) * 0.75)];
        s_p[j] = (q3 - q1) / 1.35;
      } else if (scale == "s") {
        double sum = 0.0;
        for (double v : col) sum += v;
        double mu_col = sum / n;
        double ss = 0.0;
        for (double v : col) ss += (v - mu_col) * (v - mu_col);
        s_p[j] = std::sqrt(ss / (n - 1));
      } else {
        stop("Invalid scale flag: must be 'm', 'i', or 's'");
      }
      
      if (s_p[j] == 0.0) s_p[j] = 1.0;
    } else {
      s_p[j] = 1.0;
    }
  }
  
  std::vector< std::vector<int> > ord_rank_Y(p);
  std::vector< std::vector<int> > ord_rank_mu(p);
  std::vector<int> M(p, 1);
  
  for (int j = 0; j < p; ++j) {
    if (ft[j] == 2) {
      std::vector<double> vals;
      vals.reserve(n + m);
      for (int i = 0; i < n; ++i) vals.push_back(Y(i, j));
      for (int u = 0; u < m; ++u) vals.push_back(mu(u, j));
      
      std::sort(vals.begin(), vals.end());
      vals.erase(std::unique(vals.begin(), vals.end()), vals.end());
      
      int levels = (int)vals.size();
      M[j] = levels > 1 ? levels : 1;
      
      ord_rank_Y[j].resize(n);
      ord_rank_mu[j].resize(m);
      
      for (int i = 0; i < n; ++i) {
        ord_rank_Y[j][i] =
          (int)(std::lower_bound(vals.begin(), vals.end(), Y(i, j)) - vals.begin());
      }
      for (int u = 0; u < m; ++u) {
        ord_rank_mu[j][u] =
          (int)(std::lower_bound(vals.begin(), vals.end(), mu(u, j)) - vals.begin());
      }
    }
  }
  
  NumericVector D(n * m * p);
  D.attr("dim") = IntegerVector::create(n, m, p);
  
  for (int j = 0; j < p; ++j) {
    for (int u = 0; u < m; ++u) {
      for (int i = 0; i < n; ++i) {
        double diff;
        if (ft[j] == 0) {
          diff = std::abs(Y(i, j) - mu(u, j)) / s_p[j];
        } else if (ft[j] == 1) {
          diff = (Y(i, j) != mu(u, j)) ? 1.0 : 0.0;
        } else {
          double denom = double(M[j] - 1);
          diff = (denom > 0.0) ? std::abs(ord_rank_Y[j][i] - ord_rank_mu[j][u]) / denom : 0.0;
        }
        
        D[i + n * u + n * m * j] = diff;
      }
    }
  }
  
  return D;
}


// [[Rcpp::export]]
Rcpp::List E_step(const NumericMatrix& loss_by_state,
                  const NumericMatrix& Gamma) {
  int T = loss_by_state.nrow();
  int K = loss_by_state.ncol();
  
  if (Gamma.nrow() != K || Gamma.ncol() != K)
    stop("Gamma must be K x K with K = ncol(loss_by_state)");
  
  NumericMatrix V(T, K);
  
  for (int t = 0; t < T; ++t)
    for (int j = 0; j < K; ++j)
      V(t, j) = loss_by_state(t, j);
  
  for (int t = T - 2; t >= 0; --t) {
    for (int j = 0; j < K; ++j) {
      double m = V(t + 1, 0) + Gamma(0, j);
      for (int i = 1; i < K; ++i) {
        double cand = V(t + 1, i) + Gamma(i, j);
        if (cand < m) m = cand;
      }
      V(t, j) = loss_by_state(t, j) + m;
    }
  }
  
  IntegerVector s(T);
  
  {
    double m0 = V(0, 0);
    int idx = 0;
    for (int j = 1; j < K; ++j) {
      if (V(0, j) < m0) {
        m0 = V(0, j);
        idx = j;
      }
    }
    s[0] = idx + 1;
  }
  
  for (int t = 1; t < T; ++t) {
    int prev = s[t - 1] - 1;
    double m = V(t, 0) + Gamma(prev, 0);
    int idx = 0;
    for (int j = 1; j < K; ++j) {
      double cand = V(t, j) + Gamma(prev, j);
      if (cand < m) {
        m = cand;
        idx = j;
      }
    }
    s[t] = idx + 1;
  }
  
  return Rcpp::List::create(
    Rcpp::Named("s") = s,
    Rcpp::Named("V") = V
  );
}


// [[Rcpp::export]]
Rcpp::IntegerVector pam_fixed_cpp(
    const Rcpp::NumericMatrix& D,
    const Rcpp::IntegerVector& s,
    const int K
) {
  const int TT = D.nrow();
  
  if (D.ncol() != TT)
    Rcpp::stop("D must be a square distance matrix.");
  
  if (s.size() != TT)
    Rcpp::stop("length(s) must match nrow(D).");
  
  Rcpp::IntegerVector medoids(K, NA_INTEGER);
  
  for (int k = 1; k <= K; ++k) {
    
    // collect indices of cluster k
    std::vector<int> idx;
    idx.reserve(TT);
    
    for (int t = 0; t < TT; ++t) {
      if (Rcpp::IntegerVector::is_na(s[t]))
        Rcpp::stop("s contains NA values.");
      
      if (s[t] == k)
        idx.push_back(t);
    }
    
    const int nk = idx.size();
    
    if (nk == 0) {
      medoids[k - 1] = NA_INTEGER;
      continue;
    }
    
    if (nk == 1) {
      medoids[k - 1] = idx[0] + 1; // R indexing
      continue;
    }
    
    double best_obj = R_PosInf;
    int best_j = NA_INTEGER;
    
    for (int jj = 0; jj < nk; ++jj) {
      int j = idx[jj];
      double obj = 0.0;
      
      for (int ii = 0; ii < nk; ++ii) {
        double val = D(idx[ii], j);
        if (!Rcpp::NumericVector::is_na(val))
          obj += val;
      }
      
      if (obj < best_obj) {
        best_obj = obj;
        best_j = j;
      }
    }
    
    medoids[k - 1] = best_j + 1; // back to R indexing
  }
  
  return medoids;
}

