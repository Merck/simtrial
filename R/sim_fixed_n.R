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
#' @param rho_gamma As in [tenFHcorr()]. A `tibble` with variables
#'   `rho` and `gamma`, both greater than equal to zero,
#'   to specify one Fleming-Harrington weighted logrank test per row.
#' @param seed Optional. Initial seed for simulations.
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
#' A `tibble` including columns:
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
#' @import dplyr
#' @importFrom tibble tibble
#'
#' @export
#'
#' @examples
#' library(tidyr)
#' library(dplyr)
#' library(doParallel)
#' library(tibble)
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
#'   summarise(mean(z <= qnorm(.025)))
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
#' registerDoParallel(2)
#' sim_fixed_n(n_sim = 10, seed = 2022)
#' stopImplicitCluster()
#' registerDoSEQ()
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
    rho_gamma = tibble(rho = 0, gamma = 0),
    seed = NULL) {
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

  # Check seed
  if (is.null(seed)) {
    setSeed <- FALSE
  } else {
    if (!is.numeric(seed)) {
      stop("sim_fixed_n: seed in `sim_fixed_n()` must be a number.")
    }
    setSeed <- TRUE
  }

  n_stratum <- nrow(stratum)

  # Build a function to calculate z and log-hr
  doAnalysis <- function(d, rho_gamma, n_stratum) {
    if (nrow(rho_gamma) == 1) {
      z <- tibble(z = (d %>% counting_process(arm = "experimental") %>% wlr(rho_gamma = rho_gamma))$z)
    } else {
      z <- d %>%
        counting_process(arm = "experimental") %>%
        tenFHcorr(rho_gamma = rho_gamma, corr = TRUE)
    }

    ans <- tibble(
      event = sum(d$event),
      ln_hr = ifelse(n_stratum > 1,
        survival::coxph(survival::Surv(tte, event) ~ (treatment == "experimental") + survival::strata(stratum), data = d)$coefficients,
        survival::coxph(survival::Surv(tte, event) ~ (treatment == "experimental"), data = d)$coefficients
      ) %>% as.numeric()
    )

    ans <- cbind(ans, z)
    return(ans)
  }

  # Compute minimum planned follow-up time
  minFollow <- max(0, total_duration - sum(enroll_rate$duration))

  # Put failure rates into sim_pw_surv format
  temp <- simfix2simpwsurv(fail_rate)
  fr <- temp$fail_rate
  dr <- temp$dropout_rate
  results <- NULL

  # Parallel
  `%op%` <- get_operator()
  results <- foreach::foreach(
    i = seq_len(n_sim),
    .combine = "rbind",
    .errorhandling = "pass"
  ) %op% {
    # Set random seed
    if (setSeed) {
      set.seed(seed + i - 1)
    }

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
    tedate <- sim %>% get_cut_date_by_event(target_event)

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
      d <- sim %>% cut_data_by_date(total_duration)
      r1 <- d %>% doAnalysis(rho_gamma, n_stratum)
    }

    # Targeted events cutoff
    if (tests[2]) {
      d <- sim %>% cut_data_by_date(tedate)
      r2 <- d %>% doAnalysis(rho_gamma, n_stratum)
    }

    # Minimum follow-up cutoff
    if (tests[3]) {
      d <- sim %>% cut_data_by_date(tmfdate)
      r3 <- d %>% doAnalysis(rho_gamma, n_stratum)
    }

    addit <- NULL
    # Planned duration cutoff
    if (1 %in% timing_type) {
      addit <- rbind(
        addit,
        r1 %>% mutate(
          cut = "Planned duration",
          duration = total_duration
        )
      )
    }

    # Targeted events cutoff
    if (2 %in% timing_type) {
      addit <- rbind(
        addit,
        r2 %>% mutate(
          cut = "Targeted events",
          duration = tedate
        )
      )
    }

    # Minimum follow-up duration target
    if (3 %in% timing_type) {
      addit <- rbind(
        addit,
        r3 %>% mutate(
          cut = "Minimum follow-up",
          duration = tmfdate
        )
      )
    }

    # Max of planned duration, targeted events
    if (4 %in% timing_type) {
      if (tedate > total_duration) {
        addit <- rbind(
          addit,
          r2 %>% mutate(
            cut = "Max(planned duration, event cut)",
            duration = tedate
          )
        )
      } else {
        addit <- rbind(
          addit,
          r1 %>% mutate(
            cut = "Max(planned duration, event cut)",
            duration = total_duration
          )
        )
      }
    }

    # Max of minimum follow-up, targeted events
    if (5 %in% timing_type) {
      if (tedate > tmfdate) {
        addit <- rbind(
          addit,
          r2 %>% mutate(
            cut = "Max(min follow-up, event cut)",
            duration = tedate
          )
        )
      } else {
        addit <- rbind(
          addit,
          r3 %>% mutate(
            cut = "Max(min follow-up, event cut)",
            duration = tmfdate
          )
        )
      }
    }

    addit %>% mutate(sim = i)
  }

  return(results)
}

# Get operator (parallel or serial)
get_operator <- function() {
  is_par <- foreach::getDoParWorkers() > 1
  if (is_par) {
    message("Using ", foreach::getDoParWorkers(), " cores with backend ", foreach::getDoParName())
    res <- foreach::`%dopar%`
  } else {
    res <- foreach::`%do%`
  }
  res
}
