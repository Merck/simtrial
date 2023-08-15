#  Copyright (c) 2023 Merck & Co., Inc., Rahway, NJ, USA and its affiliates.
#  All rights reserved.
#
#  This file is part of the simtrial program.
#
#  simtrial is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program.  If not, see <http://www.gnu.org/licenses/>.

#' Simulation of fixed sample size design for time-to-event endpoint
#'
#' `sim_fixed_n()` provides simulations of a single endpoint two-arm trial
#' where the enrollment, hazard ratio, and failure and dropout rates change
#' over time.
#'
#' @param n_sim Number of simulations to perform.
#' @param sample_size Total sample size per simulation.
#' @param target_event Targeted event count for analysis.
#' @param stratum A tibble with stratum specified in `stratum`,
#'   probability (incidence) of each stratum in `p`.
#' @param enroll_rate Piecewise constant enrollment rates by time period.
#'   Note that these are overall population enrollment rates and the `stratum`
#'   argument controls the random distribution between stratum.
#' @param fail_rate Piecewise constant control group failure rates, hazard ratio
#'   for experimental vs. control, and dropout rates by stratum and time period.
#' @param total_duration Total follow-up from start of enrollment to data cutoff.
#' @param block As in [sim_pw_surv()]. Vector of treatments to be included
#'   in each block.
#' @param timing_type A numeric vector determining data cutoffs used;
#'   see details. Default is to include all available cutoff methods.
#' @param rho_gamma As in [wlr()]. A `tibble` with variables
#'   `rho` and `gamma`, both greater than equal to zero,
#'   to specify one Fleming-Harrington weighted logrank test per row.
#'
#' @details
#' `timing_type` has up to 5 elements indicating different options
#' for data cutoff:
#' - `1`: Uses the planned study duration.
#' - `2`: The time the targeted event count is achieved.
#' - `3`: The planned minimum follow-up after enrollment is complete.
#' - `4`: The maximum of planned study duration and targeted event count cuts
#'   (1 and 2).
#' - `5`: The maximum of targeted event count and minimum follow-up cuts
#'   (2 and 3).
#'
#' @return
#' A data frame including columns:
#' - `event`: Event count.
#' - `ln_hr`: Log-hazard ratio.
#' - `z`: Normal test statistic; < 0 favors experimental.
#' - `cut`: Text describing cutoff used.
#' - `duration`: Duration of trial at cutoff for analysis.
#' - `sim`: Sequential simulation ID.
#'
#' One row per simulated dataset per cutoff specified in `timing_type`,
#' per test statistic specified.
#' If multiple Fleming-Harrington tests are specified in `rho_gamma`,
#' then columns `rho` and `gamma` are also included.
#'
#' @importFrom tibble tibble
#' @importFrom data.table ":=" rbindlist setDF
#' @importFrom doFuture "%dofuture%"
#' @importFrom future plan
#' @importFrom methods is
#'
#' @export
#'
#' @examples
#' library(tibble)
#' library(dplyr)
#' library(future)
#'
#' # Example 1
#' # Show output structure
#' sim_fixed_n(n_sim = 3)
#'
#' # Example 2
#' # Example with 2 tests: logrank and FH(0,1)
#' sim_fixed_n(n_sim = 1, rho_gamma = tibble(rho = 0, gamma = c(0, 1)))
#'
#' # Example 3
#' # Power by test
#' # Only use cuts for events, events + min follow-up
#' xx <- sim_fixed_n(
#'   n_sim = 100,
#'   timing_type = c(2, 5),
#'   rho_gamma = tibble(rho = 0, gamma = c(0, 1))
#' )
#' # Get power approximation for FH, data cutoff combination
#' xx %>%
#'   group_by(cut, rho, gamma) %>%
#'   summarize(mean(z <= qnorm(.025)))
#'
#' # MaxCombo power estimate for cutoff at max of targeted events, minimum follow-up
#' p <- xx %>%
#'   filter(cut != "Targeted events") %>%
#'   group_by(sim) %>%
#'   group_map(pvalue_maxcombo) %>%
#'   unlist()
#'
#' mean(p < .025)
#'
#' # MaxCombo estimate for targeted events cutoff
#' p <- xx %>%
#'   filter(cut == "Targeted events") %>%
#'   group_by(sim) %>%
#'   group_map(pvalue_maxcombo) %>%
#'   unlist()
#'
#' mean(p < .025)
#'
#' # Example 3
#' # Use two cores
#' set.seed(2023)
#' plan("multisession", workers = 2)
#' sim_fixed_n(n_sim = 10)
#' plan("sequential")
sim_fixed_n <- function(
    n_sim = 1000,
    sample_size = 500, # Sample size
    target_event = 350, # Targeted total event count
    # Multinomial probability distribution for stratum enrollment
    stratum = tibble(stratum = "All", p = 1),
    # Enrollment rates as in AHR()
    enroll_rate = tibble(duration = c(2, 2, 10), rate = c(3, 6, 9)),
    # Failure rates as in AHR()
    fail_rate = tibble(
      stratum = "All",
      duration = c(3, 100),
      fail_rate = log(2) / c(9, 18),
      hr = c(.9, .6),
      dropout_rate = rep(.001, 2)
    ),
    # Total planned trial duration; single value
    total_duration = 30,
    # Fixed block randomization specification
    block = rep(c("experimental", "control"), 2),
    # Select desired cutoffs for analysis (default is all types)
    timing_type = 1:5,
    # Default is to to logrank testing, but one or more Fleming-Harrington tests
    # can be specified
    rho_gamma = tibble(rho = 0, gamma = 0)) {
  # Check input values
  # Check input enrollment rate assumptions
  if (!("duration" %in% names(enroll_rate))) {
    stop("sim_fixed_n: enrollRates column names in `sim_fixed_n()` must contain duration.")
  }

  if (!("rate" %in% names(enroll_rate))) {
    stop("sim_fixed_n: enrollRates column names in `sim_fixed_n()` must contain rate.")
  }

  # Check input failure rate assumptions
  if (!("stratum" %in% names(fail_rate))) {
    stop("sim_fixed_n: fail_rate column names in `sim_fixed_n()` must contain stratum.")
  }

  if (!("duration" %in% names(fail_rate))) {
    stop("sim_fixed_n: fail_rate column names in `sim_fixed_n()` must contain duration.")
  }

  if (!("fail_rate" %in% names(fail_rate))) {
    stop("sim_fixed_n: fail_rate column names in `sim_fixed_n()` must contain fail_rate.")
  }

  if (!("hr" %in% names(fail_rate))) {
    stop("sim_fixed_n: fail_rate column names in `sim_fixed_n()` must contain hr")
  }

  if (!("dropout_rate" %in% names(fail_rate))) {
    stop("sim_fixed_n: fail_rate column names in `sim_fixed_n()` must contain dropout_rate")
  }

  # Check input trial duration
  if (!is.numeric(total_duration)) {
    stop("sim_fixed_n: total_duration in `sim_fixed_n()` must be a single positive number.")
  }

  if (!is.vector(total_duration)) {
    stop("sim_fixed_n: total_duration in `sim_fixed_n()` must be a single positive number.")
  }

  if (length(total_duration) != 1) {
    stop("sim_fixed_n: total_duration in `sim_fixed_n()` must be a single positive number.")
  }

  if (!min(total_duration) > 0) {
    stop("sim_fixed_n: total_duration in `sim_fixed_n()` must be a single positive number.")
  }

  # Check stratum
  stratum2 <- names(table(fail_rate$stratum))
  if (nrow(stratum) != length(stratum2)) {
    stop("sim_fixed_n: stratum in `sim_fixed_n()` must be the same in stratum and fail_rate.")
  }

  if (any(is.na(match(stratum$stratum, stratum2))) | any(is.na(match(stratum2, stratum$stratum)))) {
    stop("sim_fixed_n: stratum in `sim_fixed_n()` must be the same in stratum and fail_rate.")
  }

  # Check n_sim
  if (n_sim <= 0) {
    stop("sim_fixed_n: n_sim in `sim_fixed_n()` must be positive integer.")
  }

  if (length(n_sim) != 1) {
    stop("sim_fixed_n: n_sim in `sim_fixed_n()` must be positive integer.")
  }

  if (n_sim != ceiling(n_sim)) {
    stop("sim_fixed_n: n_sim in `sim_fixed_n()` must be positive integer.")
  }

  # Check target_event
  if (target_event <= 0) {
    stop("sim_fixed_n: target_event in `sim_fixed_n()` must be positive.")
  }

  if (length(target_event) != 1) {
    stop(("sim_fixed_n: target_event in `sim_fixed_n()` must be positive."))
  }

  # Check sample_size
  if (sample_size <= 0) {
    stop("sim_fixed_n: sample_size in `sim_fixed_n()` must be positive.")
  }

  if (length(sample_size) != 1) {
    stop("sim_fixed_n: sample_size in `sim_fixed_n()` must be positive.")
  }

  n_stratum <- nrow(stratum)

  # Compute minimum planned follow-up time
  minFollow <- max(0, total_duration - sum(enroll_rate$duration))

  # Put failure rates into sim_pw_surv format
  temp <- simfix2simpwsurv(fail_rate)
  fr <- temp$fail_rate
  dr <- temp$dropout_rate
  results <- NULL

  # message for backends
  if (!is(plan(), "sequential")) {
    # future backend
    message("Using ", nbrOfWorkers(), " cores with backend ", attr(plan("list")[[1]], "class")[2])
  } else if (foreach::getDoParWorkers() > 1) {
    message("Using ", foreach::getDoParWorkers(), " cores with backend ", foreach::getDoParName())
    message("Warning: ")
    message("doFuture may exhibit suboptimal performance when using a doParallel backend.")
  } else {
    message("Backend uses sequential processing.")
  }

  results <- foreach::foreach(
    i = seq_len(n_sim),
    .combine = "rbind",
    .errorhandling = "pass",
    .options.future = list(seed = TRUE)
  ) %dofuture% {
    # Generate piecewise data
    sim <- sim_pw_surv(
      n = sample_size,
      stratum = stratum,
      enroll_rate = enroll_rate,
      fail_rate = fr,
      dropout_rate = dr,
      block = block
    )

    # Study date that targeted event rate achieved
    tedate <- get_cut_date_by_event(sim, target_event)

    # Study data that targeted minimum follow-up achieved
    tmfdate <- max(sim$enroll_time) + minFollow

    # Compute tests for all specified cutoff options
    r1 <- NULL
    r2 <- NULL
    r3 <- NULL

    tests <- rep(FALSE, 3)

    ## Total planned trial duration or max of this and targeted events
    if (1 %in% timing_type) {
      tests[1] <- TRUE # Planned duration cutoff
    }

    if (2 %in% timing_type) {
      tests[2] <- TRUE # Targeted events cutoff
    }

    if (3 %in% timing_type) {
      tests[3] <- TRUE # Minimum follow-up duration target
    }

    if (4 %in% timing_type) { # Max of planned duration, targeted events
      if (tedate > total_duration) {
        tests[2] <- TRUE
      } else {
        tests[1] <- TRUE
      }
    }

    if (5 %in% timing_type) { # Max of minimum follow-up, targeted events
      if (tedate > tmfdate) {
        tests[2] <- TRUE
      } else {
        tests[3] <- TRUE
      }
    }

    # Total duration cutoff
    if (tests[1]) {
      d <- cut_data_by_date(sim, total_duration)
      r1 <- doAnalysis(d, rho_gamma, n_stratum)
    }

    # Targeted events cutoff
    if (tests[2]) {
      d <- cut_data_by_date(sim, tedate)
      r2 <- doAnalysis(d, rho_gamma, n_stratum)
    }

    # Minimum follow-up cutoff
    if (tests[3]) {
      d <- cut_data_by_date(sim, tmfdate)
      r3 <- doAnalysis(d, rho_gamma, n_stratum)
    }

    addit <- list()
    # Planned duration cutoff
    if (1 %in% timing_type) {
      rtemp <- cbind(r1, cut = "Planned duration", duration = total_duration)
      addit <- c(addit, list(rtemp))
    }

    # Targeted events cutoff
    if (2 %in% timing_type) {
      rtemp <- cbind(r2, cut = "Targeted events", duration = tedate)
      addit <- c(addit, list(rtemp))
    }

    # Minimum follow-up duration target
    if (3 %in% timing_type) {
      rtemp <- cbind(r3, cut = "Minimum follow-up", duration = tmfdate)
      addit <- c(addit, list(rtemp))
    }

    # Max of planned duration, targeted events
    if (4 %in% timing_type) {
      if (tedate > total_duration) {
        rtemp <- cbind(r2, cut = "Max(planned duration, event cut)", duration = tedate)
        addit <- c(addit, list(rtemp))
      } else {
        rtemp <- cbind(r1, cut = "Max(planned duration, event cut)", duration = total_duration)
        addit <- c(addit, list(rtemp))
      }
    }

    # Max of minimum follow-up, targeted events
    if (5 %in% timing_type) {
      if (tedate > tmfdate) {
        rtemp <- cbind(r2, cut = "Max(min follow-up, event cut)", duration = tedate)
        addit <- c(addit, list(rtemp))
      } else {
        rtemp <- cbind(r3, cut = "Max(min follow-up, event cut)", duration = tmfdate)
        addit <- c(addit, list(rtemp))
      }
    }

    results_sim <- rbindlist(addit)
    results_sim[, sim := i]
    setDF(results_sim)
    return(results_sim)
  }

  return(results)
}

# Build a function to calculate z and log-hr
doAnalysis <- function(d, rho_gamma, n_stratum) {
  if (nrow(rho_gamma) == 1) {
    tmp <- counting_process(d, arm = "experimental")
    tmp <- wlr(tmp, rho_gamma = rho_gamma)
    z <- data.frame(z = tmp$z)
  } else {
    tmp <- counting_process(d, arm = "experimental")
    z <- wlr(tmp, rho_gamma = rho_gamma, return_corr = TRUE)
  }

  event <- sum(d$event)
  if (n_stratum > 1) {
    ln_hr <- survival::coxph(survival::Surv(tte, event) ~ (treatment == "experimental") + survival::strata(stratum), data = d)$coefficients
  } else {
    ln_hr <- survival::coxph(survival::Surv(tte, event) ~ (treatment == "experimental"), data = d)$coefficients
  }
  ln_hr <- as.numeric(ln_hr)

  ans <- cbind(event, ln_hr, z)
  return(ans)
}
