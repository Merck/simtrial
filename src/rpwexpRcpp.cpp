#include <Rcpp.h>
#include <random>
using namespace Rcpp;

#define SIMTRIAL_DEBUG

//' @param n Number of observations to be generated.
//' @param failRates A dataframe containing \code{duration} and \code{rate} variables.
//' @export
// [[Rcpp::export]]
NumericVector rpwexpRcpp(int n,
                      DataFrame failRates) {
  NumericVector duration = failRates["duration"];
  NumericVector rate = failRates["rate"];
  int n_rates = duration.size();

  // Initialize failure times to Inf
  NumericVector times(n, R_PosInf);

#ifdef SIMTRIAL_DEBUG
  // Set seed for comparison
  Environment base_env("package:base");
  Function set_seed_r = base_env["set.seed"];
  set_seed_r(2022);
#endif

  if (n_rates == 1) {
    if (rate[0] != 0) {
      times = rexp(n, rate[0]); // generate exponential failure time if non-0 failure rate
    }
  } else {
    double starttime = 0; // start of first failure rate interval
    LogicalVector indx(n, 1); // index for event times not yet reached
    int nindx = n; // number of event times left to generate

    for (int i = 0; i < n_rates; i++) {
      if (is_false(any(indx))) break; // stop if finished
      if (rate[i] == 0) {
        NumericVector temp(nindx, R_PosInf); // set failure time to Inf for inveral i if 0 fail rate
        times[indx] = temp;
      } else {
        NumericVector temp = rexp(nindx, rate[i]); // generate exponential failure time for interval i if non-0 faiurel rate
        int p = 0;
        for (int j = 0; j < n; j++) {
          if (indx[j]) times[j] = temp[p++] + starttime;
        }
      }

      if (i < n_rates) { // skip this for last interval as all remaining times are generated there
        starttime += duration[i]; // update start time for next interval
        for (int j = 0; j < n; j++) {
          if (indx[j] == 1 && times[j] <= starttime) {
            indx[j] = 0; // update index of event times not yet reached
            nindx--; // update number of event times left to generate
          }
        }
      }

    }
  }

  return times;
}
