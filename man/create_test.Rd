% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sim_gs_n.R
\name{create_test}
\alias{create_test}
\title{Create a cutting test function}
\usage{
create_test(test, ...)
}
\arguments{
\item{test}{A test function such as \code{\link[=wlr]{wlr()}}, \code{\link[=maxcombo]{maxcombo()}}, or \code{\link[=rmst]{rmst()}}}

\item{...}{Arguments passed to the cutting test function}
}
\value{
A function that accepts a data frame of simulated trial data and
returns a test result
}
\description{
Create a cutting test function for use with \code{\link[=sim_gs_n]{sim_gs_n()}}
}
\examples{
# Simulate trial data
trial_data <- sim_pw_surv()

# Cut after 150 events
trial_data_cut <- cut_data_by_event(trial_data, 150)

# Create a cutting test function that can be used by sim_gs_n()
regular_logrank_test <- create_test(wlr, weight = fh(rho = 0, gamma = 0))

# Test the cutting
regular_logrank_test(trial_data_cut)

# The results are the same as directly calling the function
stopifnot(all.equal(
  regular_logrank_test(trial_data_cut),
  wlr(trial_data_cut, weight = fh(rho = 0, gamma = 0))
))
}
\seealso{
\code{\link[=sim_gs_n]{sim_gs_n()}}, \code{\link[=create_cut]{create_cut()}}
}
