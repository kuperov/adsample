// Adaptive rejection sampling algorithm for R based on Gilks & Wild (1992).
// Allows efficient sampling from univariate densities that are log concave.
//
// Alex Cooper <alexander.cooper@monash.edu>
// August 2017

#include <cmath>
#include <cstdlib>  // for abs()
#include <vector>
#include <utility>

using namespace std;

#include <Rcpp.h>
using namespace Rcpp;

// (point value, gradient value)
typedef pair<double, double> pointgrad;

// callback to evaluate the (possibly expensive) log density function,
// which returns the ordinate and the gradient evaluated at the abscissa x.
// Note that while the gradient need not be normalized (it can differ from
// the true gradient up to a constant multiple), the gradient and ordinate
// must share the *same* constant.
typedef std::function<pointgrad(double x)> log_dens_callback;

// container for the state of the algorithm, for returning state back to
// R for debugging
typedef struct {
  int k;
  vector<double> T;
  vector<double> H;
  vector<double> Hprime;
  vector<double> Z;
  vector<double> Q;
} algo_state;

// Draw n variates from the unnormalized density h using adaptive sampling.
//
// IMPORTANT: the density MUST be univariate and log concave in this parameter.
//
// params:
//    n - number of variates to draw
//    state - optional pointer to a map for storing algorithm state. Leave as
//            NULL to discard state at end of algorithm
//    maxiter - maximum number of iterations, zero for no limit
vector<double> adsample(int n,
                        log_dens_callback h,
                        std::vector<double> initAbscissae,
                        double minD,
                        double maxD,
                        algo_state* state = NULL,
                        long maxiter = 0) {
  vector<double> variates(n);

  int k = initAbscissae.size();

  vector<double> T(initAbscissae); // [x_1, x_2, ..., x_k]
  vector<double> H(k);             // [h(x_1), h(x_2), ..., h(x_k)]
  vector<double> Hprime(k);        // [h'(x_1), h'(x_2), ..., h'(x_k)]
  for (int i=0; i<k; i++) {
    auto pgrad = h(T[i]);
    H[i] = pgrad.first;
    Hprime[i] = pgrad.second;
  }

  // Compute the value of z_i. This is equation (1) in the paper,
  // and unlike the other variables in this program the subscript used
  // in the paper lines up exactly with the i parameter used here.
  auto z = [&T, &H, &Hprime](int i) {
    return (H[i] - H[i-1] - T[i]*Hprime[i] + T[i-1]*Hprime[i-1])/
      (Hprime[i-1] - Hprime[i]);
  };

  // partition of piecewise envelope function
  vector<double> Z(k + 1);
  Z[0] = minD;
  Z[k] = maxD;
  for (int i=1; i < k; i++) {
    Z[i] = z(i);
  }

  // Compute the area under the unnormalized envelope curve Q*s_k(x) between
  // z_i and z_{i+1}, ie the (i+1)th piecewise segment of Q*s_k(k).
  //
  // Note this is not the normalized s_k(x) given in equation (3), but instead
  //     Q*s_k(x) = exp(u(x)) = exp(h(x_1) + (x - x_1)*h'(x)),
  // ie Q has the value of the integral over D following the / in equation (3).
  //
  auto env_area_pt = [&T, &H, &Hprime, &Z](int i) {
    return exp(H[i] - T[i]*Hprime[i])*
      (exp(Z[i+1]*Hprime[i]) - exp(Z[i]*Hprime[i]))/Hprime[i];
  };

  // Qtot and Q are area under s_k(x) partitioned by Zs
  vector<double> Q(k);
  double Qtot = 0.;
  for (int i=0; i<k; i++) {
    double area = env_area_pt(i);
    Q[i] = area;
    Qtot += area;
  }

  int nsampled = 0;
  for (long iter = 0; nsampled < n && (maxiter == 0 || iter >= maxiter);
       iter++) {
    // step 1: draw from s_k(x)
    double w = as<double>(runif(1));
    // renormalize this draw to the measure of envelope s_k(x)
    double wq = Qtot * w;
    int i = 0;  // index for piecewise chunk of s_k(x)
    double q_partial = 0; // area to the left of ith chunk
    while (q_partial + Q[i] < wq)
      q_partial += Q[i++];
    // xstar is the x ordinate corresponding to w
    double xstar = log((wq - q_partial) * Hprime[i] * exp(T[i]*Hprime[i] - H[i])
                       + exp(Z[i] * Hprime[i])) / Hprime[i];

    // step 2: squeezing test (sec 2.2.2) using cheap upper and lower hulls
    w = as<double>(runif(1));
    double u_xstar = H[i] + (xstar - T[i])*Hprime[i];  // equation (2)
    double l_xstar = ((i > 0) && (i < k-1)) ?
    ((T[i+1] - xstar)*H[i] + (xstar - T[i])*H[i+1])/(T[i+1] - T[i])  // eq (4)
      :
      -std::numeric_limits<double>::infinity();   // outside support of l(x)
    if (w <= exp(l_xstar - u_xstar)) {
      variates[nsampled++] = xstar;
      continue;
    }

    // squeezing test has failed, so we now pay the cost of evaluating h(x)
    // and h'(x), and will add an additional abscissa to T below
    auto h_xstar_pair = h(xstar);
    double h_xstar = h_xstar_pair.first;
    double h_xstar_prime = h_xstar_pair.second;
    // rejection test (sec 2.2.2 contd)
    if (w <= exp(h_xstar - u_xstar)) {
      variates[nsampled++] = xstar;
      // no continue here since we still need to save h(xstar) and h'(xstar)
    }

    // step 3: save h(xstar) and h'(xstar) and recompute the others
    k += 1;
    if (xstar < T[i]) {
      // xstar inserted to LEFT of x_j
      int j = i;

      T.insert(T.begin()+j, xstar);
      H.insert(H.begin()+j, h_xstar);
      Hprime.insert(Hprime.begin()+j, h_xstar_prime);
      Z.insert(Z.begin()+(j+1), z(j+1));
      if (j != 0) { // LH boundary is never updated
        Z[j] = z(j);
      }
      Qtot -= Q[i];
      Q.insert(Q.begin()+j, env_area_pt(j));
      // we redraw the boundary of region to right
      Q[j + 1] = env_area_pt(j+1);
      Qtot += Q[j + 1] + Q[j];
      if (j != 0) {
        Qtot -= Q[j-1];
        // we redraw the boundary of region to left
        Q[j-1] = env_area_pt(j-1);
        Qtot += Q[j-1];
      }
    } else {
      // xstar inserted to RIGHT of x_j
      int j = i + 1;
      /* printf("RHS j = %d\n", j); */
      T.insert(T.begin()+j, xstar);
      H.insert(H.begin()+j, h_xstar);
      Hprime.insert(Hprime.begin()+j, h_xstar_prime);
      Z.insert(Z.begin()+j, z(j));
      if (j + 1 < k) {  // RH boundary is never updated
        Z[j+1] = z(j+1);
      }
      // Z[j] = z(j, T, H, Hprime);
      Qtot -= Q[i];
      Q.insert(Q.begin()+j, env_area_pt(j));
      // we redraw the boundary of region to left
      Q[j-1] = env_area_pt(j-1);
      Qtot += Q[j] + Q[j-1];
      if (j < k-1) {
        Qtot -= Q[j+1];
        // we redraw the boundary of region to right
        Q[j+1] = env_area_pt(j+1);
        Qtot += Q[j+1];
      }
    }
  }

  if (state != NULL) {
    state->k = k;
    state->T = T;
    state->H = H;
    state->Hprime = Hprime;
    state->Z = Z;
    state->Q = Q;
  }

  return variates;
}

//' Rcpp wrapper for adsamp algorithm
//'
//' @param n number of variates to draw
//' @param log_dens log density function, a function with a single parameter
//' @param initialPoints at least 2 points in support to seed algorithm
//' @param minRange lower bound of support
//' @param maxRange upper bound of support
//'
// [[Rcpp::export]]
NumericVector raw_ad_sample(int n,
                            Rcpp::Function log_dens,
                            NumericVector initialPoints,
                            Rcpp::NumericVector minRange,
                            Rcpp::NumericVector maxRange) {

  // lambda for getting log density and its slope
  log_dens_callback h = [log_dens](double x) {
    NumericVector pointslope = log_dens(x);
    double y = pointslope[0];
    double yprime = pointslope[1];
    return std::pair<double, double>(y, yprime);
  };
  // initial abscissae
  std::vector<double> initT(initialPoints.begin(), initialPoints.end());
  // boundaries of the support of h
  double minD = Rcpp::as<double>(minRange);
  double maxD = Rcpp::as<double>(maxRange);

  auto samp = adsample(n, h, initT, minD, maxD, NULL, 0);
  return NumericVector(samp.begin(), samp.end());
}

//' Duplicate wrapper for adsamp algorithm for debugging. This version
//' returns the algorithm state as well as the variates.
//'
//' @param n number of variates to draw
//' @param log_dens log density function, a function with a single parameter
//' @param initialPoints at least 2 points in support to seed algorithm
//' @param minRange lower bound of support
//' @param maxRange upper bound of support
//' @param maxiter maximum number of iterations; zero if no limit
//'
// [[Rcpp::export]]
List raw_ad_sample_debug(int n,
                         Rcpp::Function log_dens,
                         NumericVector initialPoints,
                         Rcpp::NumericVector minRange,
                         Rcpp::NumericVector maxRange,
                         long maxiter) {

  // lambda for getting log density and its slope
  log_dens_callback h = [log_dens](double x) {
    NumericVector pointslope = log_dens(x);
    double y = pointslope[0];
    double yprime = pointslope[1];
    return std::pair<double, double>(y, yprime);
  };
  // initial abscissae
  std::vector<double> initT(initialPoints.begin(), initialPoints.end());
  // boundaries of the support of h
  double minD = Rcpp::as<double>(minRange);
  double maxD = Rcpp::as<double>(maxRange);

  algo_state state;
  auto samp = adsample(n, h, initT, minD, maxD, &state, maxiter);
  auto variates = NumericVector(samp.begin(), samp.end());

  List results;
  results["samples"] = variates;
  results["T"] = NumericVector(state.T.begin(), state.T.end());
  results["k"] = state.k;
  results["Z"] = state.Z;
  results["H"] = state.H;
  results["Hprime"] = state.Hprime;
  results["Q"] = state.Q;
  return results;
}
