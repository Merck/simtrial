#include <Rcpp.h>
#include <random>
using namespace Rcpp;

//' The piecewise exponential distribution in C++
//'
//' @param n Number of observations to be generated.
//' @param fail_rate A data frame containing `duration` and `rate` variables.
//'
//' @export
// [[Rcpp::export]]
NumericVector rpwexpRcpp(int n, DataFrame fail_rate)
{
  NumericVector duration = fail_rate["duration"];
  NumericVector rate = fail_rate["rate"];
  int n_rates = duration.size();

  // Initialize failure times to Inf
  NumericVector times(n, R_PosInf);

  if (n_rates == 1)
  {
    if (rate[0] != 0)
    {
      times = rexp(n, rate[0]); // Generate exponential failure time if non-zero failure rate
    }
  }
  else
  {
    double starttime = 0;     // Start of first failure rate interval
    LogicalVector indx(n, 1); // Index for event times not yet reached
    int nindx = n;            // Number of event times left to generate

    for (int i = 0; i < n_rates; i++)
    {
      if (is_false(any(indx)))
        break; // stop if finished
      if (rate[i] == 0)
      {
        NumericVector temp(nindx, R_PosInf); // Set failure time to Inf for interval i if 0 fail rate
        times[indx] = temp;
      }
      else
      {
        NumericVector temp = rexp(nindx, rate[i]); // Generate exponential failure time for interval i if non-0 failure rate
        int p = 0;
        for (int j = 0; j < n; j++)
        {
          if (indx[j])
            times[j] = temp[p++] + starttime;
        }
      }

      if (i < n_rates)
      {                           // Skip this for last interval as all remaining times are generated there
        starttime += duration[i]; // Update start time for next interval
        for (int j = 0; j < n; j++)
        {
          if (indx[j] == 1 && times[j] <= starttime)
          {
            indx[j] = 0; // Update index of event times not yet reached
            nindx--;     // Update number of event times left to generate
          }
        }
      }
    }
  }

  return times;
}
