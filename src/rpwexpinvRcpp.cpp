#include <Rcpp.h>
#include <random>
#include <algorithm>
using namespace Rcpp;

// #define SIMTRIAL_DEBUG

//' @param n Number of observations to be generated.
//' @param failRates A dataframe containing \code{duration} and \code{rate} variables.
//' @export
// [[Rcpp::export]]
NumericVector rpwexpinvRcpp(int n,
                            DataFrame failRates) {

  NumericVector duration = failRates["duration"];
  NumericVector rate = failRates["rate"];
  int n_rates = duration.size();

#ifdef SIMTRIAL_DEBUG
  // Set seed for comparison
  Environment base_env("package:base");
  Function set_seed_r = base_env["set.seed"];
  set_seed_r(2022);
#endif

  // Generate cumulative hazard for each observation
  NumericVector times = -log(runif(n));

  // Get number of piecewise rates
  NumericVector cumTime(n_rates);
  NumericVector cumHaz(n_rates);
  for (int i = 1; i < n_rates; i++) {
    cumTime[i] = cumTime[i-1] + duration[i-1];
    cumHaz[i] = cumHaz[i-1] + duration[i-1] * rate[i-1];
  }

  NumericVector::iterator pos;
  int j;
  for (int i = 0; i < n; i++) {
    pos = std::upper_bound(cumHaz.begin(), cumHaz.end(), times[i]);
    j = pos - cumHaz.begin() - 1;
    times[i] = cumTime[j] + (times[i] - cumHaz[j]) / rate[j];
  }

  return times;
}
