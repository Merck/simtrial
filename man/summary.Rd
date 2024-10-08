% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/summary.R
\name{summary.simtrial_gs_wlr}
\alias{summary.simtrial_gs_wlr}
\title{Summary of group sequential simulations.}
\usage{
\method{summary}{simtrial_gs_wlr}(object, design = NULL, bound = NULL, ...)
}
\arguments{
\item{object}{Simulation results generated by \code{\link[=sim_gs_n]{sim_gs_n()}}}

\item{design}{Asymptotic design generated by \code{\link[gsDesign2:gs_design_ahr]{gsDesign2::gs_design_ahr()}},
\code{\link[gsDesign2:gs_power_ahr]{gsDesign2::gs_power_ahr()}}, \code{\link[gsDesign2:gs_design_wlr]{gsDesign2::gs_design_wlr()}}, or \link[gsDesign2:gs_power_wlr]{gsDesign2::gs_power_wlr}.}

\item{bound}{The boundaries.}

\item{...}{Additional parameters (not used).}
}
\value{
A data frame
}
\description{
Summary of group sequential simulations.
}
\examples{
library(gsDesign2)

# Parameters for enrollment
enroll_rampup_duration <- 4 # Duration for enrollment ramp up
enroll_duration <- 16 # Total enrollment duration
enroll_rate <- define_enroll_rate(
  duration = c(
    enroll_rampup_duration, enroll_duration - enroll_rampup_duration),
 rate = c(10, 30))

# Parameters for treatment effect
delay_effect_duration <- 3 # Delay treatment effect in months
median_ctrl <- 9 # Survival median of the control arm
median_exp <- c(9, 14) # Survival median of the experimental arm
dropout_rate <- 0.001
fail_rate <- define_fail_rate(
  duration = c(delay_effect_duration, 100),
  fail_rate = log(2) / median_ctrl,
  hr = median_ctrl / median_exp,
  dropout_rate = dropout_rate)

# Other related parameters
alpha <- 0.025 # Type I error
beta <- 0.1 # Type II error
ratio <- 1 # Randomization ratio (experimental:control)

# Build a one-sided group sequential design
design <- gs_design_ahr(
  enroll_rate = enroll_rate, fail_rate = fail_rate,
  ratio = ratio, alpha = alpha, beta = beta,
  analysis_time = c(12, 24, 36),
  upper = gs_spending_bound,
  upar = list(sf = gsDesign::sfLDOF, total_spend = alpha),
  lower = gs_b,
  lpar = rep(-Inf, 3))

# Define cuttings of 2 IAs and 1 FA
ia1_cut <- create_cut(target_event_overall = ceiling(design$analysis$event[1]))
ia2_cut <- create_cut(target_event_overall = ceiling(design$analysis$event[2]))
fa_cut <- create_cut(target_event_overall = ceiling(design$analysis$event[3]))

# Run simulations
simulation <- sim_gs_n(
  n_sim = 3,
  sample_size = ceiling(design$analysis$n[3]),
  enroll_rate = design$enroll_rate,
  fail_rate = design$fail_rate,
  test = wlr,
  cut = list(ia1 = ia1_cut, ia2 = ia2_cut, fa = fa_cut),
  weight = fh(rho = 0, gamma = 0.5))

# Summarize simulations
bound <- gsDesign::gsDesign(k = 3, test.type = 1, sfu = gsDesign::sfLDOF)$upper$bound
simulation |> summary(bound = bound)

# Summarize simulation and compare with the planned design
simulation |> summary(design = design)
}
