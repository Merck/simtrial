% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pvalue_maxcombo.R
\name{pvalue_maxcombo}
\alias{pvalue_maxcombo}
\title{MaxCombo p-value}
\usage{
pvalue_maxcombo(
  z,
  algorithm = mvtnorm::GenzBretz(maxpts = 50000, abseps = 1e-05)
)
}
\arguments{
\item{z}{A dataset output from \code{\link[=fh_weight]{fh_weight()}}; see examples.}

\item{algorithm}{This is passed directly to the \code{algorithm} argument
in \code{\link[mvtnorm:pmvnorm]{mvtnorm::pmvnorm()}}.}
}
\value{
A numeric p-value.
}
\description{
Computes p-values for the MaxCombo test based on output from \code{\link[=fh_weight]{fh_weight()}}.
This is still in an experimental stage and is intended for use with
the \code{\link[=sim_fixed_n]{sim_fixed_n()}} trial simulation routine.
However, it can also be used to analyze clinical trial data such as
that provided in the ADaM ADTTE format.
}
\examples{
library(dplyr)

# Example 1
x <- sim_fixed_n(
  n_sim = 1,
  timing_type = 5,
  rho_gamma = data.frame(
    rho = c(0, 0, 1),
    gamma = c(0, 1, 1)
  )
)
head(x)
pvalue_maxcombo(x)

# Example 2
# Only use cuts for events, events + min follow-up
xx <- sim_fixed_n(
  n_sim = 100,
  timing_type = 5,
  rho_gamma = data.frame(
    rho = c(0, 0, 1),
    gamma = c(0, 1, 1)
  )
)
head(xx)

# MaxCombo power estimate for cutoff at max of targeted events, minimum follow-up
p <- xx |>
  group_by(sim) |>
  group_map(~ pvalue_maxcombo(.x)) |>
  unlist()
mean(p < .025)
}
