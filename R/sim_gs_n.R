#  Copyright (c) 2024 Merck & Co., Inc., Rahway, NJ, USA and its affiliates.
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

#' Simulate group sequential designs with fixed sample size
#'
#' This function uses the option "stop" for the error-handling behavior of the
#' foreach loop. This will cause the entire function to stop when errors are
#' encountered and return the first error encountered instead of returning
#' errors for each individual simulation.
#'
#' WARNING: This experimental function is a work-in-progress. The function
#' arguments will change as we add additional features.
#'
#' @inheritParams sim_fixed_n
#' @param test One or more test functions such as [wlr()], [rmst()], or
#'   [milestone()] ([maxcombo()] can only be applied by itself). If a single
#'   test function is provided, it will be applied at each cut. Alternatively a
#'   list of functions created by [create_test()]. The list form is experimental
#'   and currently limited. It only accepts one test per cutting (in the future
#'   multiple tests may be accepted), and all the tests must consistently return
#'   the same exact results (again this may be more flexible in the future).
#'   Importantly, note that the simulated data set is always passed as the first
#'   positional argument to each test function provided.
#' @param cut A list of cutting functions created by [create_cut()], see
#'   examples.
#' @param original_design A design object from the gsDesign2 package, which is required when users
#' want to calculate updated bounds. The default is NULL leaving the updated bounds uncalculated.
#' @param ia_alpha_spending Spend alpha at interim analysis based on
#' - `"min_planned_actual"`: the minimal of planned and actual alpha spending.
#' - `"actual"`: the actual alpha spending.
#' @param fa_alpha_spending If targeted final event count is not achieved (under-running at final analysis),
#' specify how to do final spending. Generally, this should be specified in analysis plan.
#' - `"info_frac"` = spend final alpha according to final information fraction
#' - `"full_alpha"` = spend full alpha at final analysis.
#' @param ... Arguments passed to the test function(s) provided by the argument
#'   `test`.
#'
#' @return A data frame summarizing the simulation ID, analysis date,
#'   z statistics or p-values.
#'
#' @importFrom data.table rbindlist setDF
#'
#' @export
#'
#' @examplesIf requireNamespace("gsDesign2", quietly = TRUE)
#' library(gsDesign2)
#'
#' # Parameters for enrollment
#' enroll_rampup_duration <- 4 # Duration for enrollment ramp up
#' enroll_duration <- 16 # Total enrollment duration
#' enroll_rate <- define_enroll_rate(
#'   duration = c(
#'     enroll_rampup_duration,
#'     enroll_duration - enroll_rampup_duration
#'   ),
#'   rate = c(10, 30)
#' )
#'
#' # Parameters for treatment effect
#' delay_effect_duration <- 3 # Delay treatment effect in months
#' median_ctrl <- 9 # Survival median of the control arm
#' median_exp <- c(9, 14) # Survival median of the experimental arm
#' dropout_rate <- 0.001
#' fail_rate <- define_fail_rate(
#'   duration = c(delay_effect_duration, 100),
#'   fail_rate = log(2) / median_ctrl,
#'   hr = median_ctrl / median_exp,
#'   dropout_rate = dropout_rate
#' )
#'
#' # Other related parameters
#' alpha <- 0.025 # Type I error
#' beta <- 0.1 # Type II error
#' ratio <- 1 # Randomization ratio (experimental:control)
#'
#' # Define cuttings of 2 IAs and 1 FA
#' # IA1
#' # The 1st interim analysis will occur at the later of the following 3 conditions:
#' # - At least 20 months have passed since the start of the study.
#' # - At least 100 events have occurred.
#' # - At least 20 months have elapsed after enrolling 200/400 subjects, with a
#' #   minimum of 20 months follow-up.
#' # However, if events accumulation is slow, we will wait for a maximum of 24 months.
#' ia1_cut <- create_cut(
#'   planned_calendar_time = 20,
#'   target_event_overall = 100,
#'   max_extension_for_target_event = 24,
#'   min_n_overall = 200,
#'   min_followup = 20
#' )
#'
#' # IA2
#' # The 2nd interim analysis will occur at the later of the following 3 conditions:
#' # - At least 32 months have passed since the start of the study.
#' # - At least 250 events have occurred.
#' # - At least 10 months after IA1.
#' # However, if events accumulation is slow, we will wait for a maximum of 34 months.
#' ia2_cut <- create_cut(
#'   planned_calendar_time = 32,
#'   target_event_overall = 200,
#'   max_extension_for_target_event = 34,
#'   min_time_after_previous_analysis = 10
#' )
#'
#' # FA
#' # The final analysis will occur at the later of the following 2 conditions:
#' # - At least 45 months have passed since the start of the study.
#' # - At least 300 events have occurred.
#' fa_cut <- create_cut(
#'   planned_calendar_time = 45,
#'   target_event_overall = 350
#' )
#'
#' # Example 1: regular logrank test at all 3 analyses
#' sim_gs_n(
#'   n_sim = 3,
#'   sample_size = 400,
#'   enroll_rate = enroll_rate,
#'   fail_rate = fail_rate,
#'   test = wlr,
#'   cut = list(ia1 = ia1_cut, ia2 = ia2_cut, fa = fa_cut),
#'   weight = fh(rho = 0, gamma = 0)
#' )
#'
#' # Example 2: weighted logrank test by FH(0, 0.5) at all 3 analyses
#' sim_gs_n(
#'   n_sim = 3,
#'   sample_size = 400,
#'   enroll_rate = enroll_rate,
#'   fail_rate = fail_rate,
#'   test = wlr,
#'   cut = list(ia1 = ia1_cut, ia2 = ia2_cut, fa = fa_cut),
#'   weight = fh(rho = 0, gamma = 0.5)
#' )
#'
#' # Example 3: weighted logrank test by MB(3) at all 3 analyses
#' sim_gs_n(
#'   n_sim = 3,
#'   sample_size = 400,
#'   enroll_rate = enroll_rate,
#'   fail_rate = fail_rate,
#'   test = wlr,
#'   cut = list(ia1 = ia1_cut, ia2 = ia2_cut, fa = fa_cut),
#'   weight = mb(delay = 3)
#' )
#'
#' # Example 4: weighted logrank test by early zero (6) at all 3 analyses
#' sim_gs_n(
#'   n_sim = 3,
#'   sample_size = 400,
#'   enroll_rate = enroll_rate,
#'   fail_rate = fail_rate,
#'   test = wlr,
#'   cut = list(ia1 = ia1_cut, ia2 = ia2_cut, fa = fa_cut),
#'   weight = early_zero(6)
#' )
#'
#' # Example 5: RMST at all 3 analyses
#' sim_gs_n(
#'   n_sim = 3,
#'   sample_size = 400,
#'   enroll_rate = enroll_rate,
#'   fail_rate = fail_rate,
#'   test = rmst,
#'   cut = list(ia1 = ia1_cut, ia2 = ia2_cut, fa = fa_cut),
#'   tau = 20
#' )
#'
#' # Example 6: Milestone at all 3 analyses
#' sim_gs_n(
#'   n_sim = 3,
#'   sample_size = 400,
#'   enroll_rate = enroll_rate,
#'   fail_rate = fail_rate,
#'   test = milestone,
#'   cut = list(ia1 = ia1_cut, ia2 = ia2_cut, fa = fa_cut),
#'   ms_time = 10
#' )
#'
#' # Warning: this example will be executable when we add info info0 to the milestone test
#' # Example 7: WLR with fh(0, 0.5) test at IA1,
#' # WLR with mb(6, Inf) at IA2, and milestone test at FA
#' ia1_test <- create_test(wlr, weight = fh(rho = 0, gamma = 0.5))
#' ia2_test <- create_test(wlr, weight = mb(delay = 6, w_max = Inf))
#' fa_test <- create_test(milestone, ms_time = 10)
#' \dontrun{
#' sim_gs_n(
#'   n_sim = 3,
#'   sample_size = 400,
#'   enroll_rate = enroll_rate,
#'   fail_rate = fail_rate,
#'   test = list(ia1 = ia1_test, ia2 = ia2_test, fa = fa_test),
#'   cut = list(ia1 = ia1_cut, ia2 = ia2_cut, fa = fa_cut)
#' )
#' }
#'
#' # WARNING: Multiple tests per cut will be enabled in a future version.
#' #          Currently does not work.
#' # Example 8: At IA1, we conduct 3 tests, LR, WLR with fh(0, 0.5), and RMST test.
#' # At IA2, we conduct 2 tests, LR and WLR with early zero (6).
#' # At FA, we conduct 2 tests, LR and milestone test.
#' ia1_test <- list(
#'   test1 = create_test(wlr, weight = fh(rho = 0, gamma = 0)),
#'   test2 = create_test(wlr, weight = fh(rho = 0, gamma = 0.5)),
#'   test3 = create_test(rmst, tau = 20)
#' )
#' ia2_test <- list(
#'   test1 = create_test(wlr, weight = fh(rho = 0, gamma = 0)),
#'   test2 = create_test(wlr, weight = early_zero(6))
#' )
#' fa_test <- list(
#'   test1 = create_test(wlr, weight = fh(rho = 0, gamma = 0)),
#'   test3 = create_test(milestone, ms_time = 20)
#' )
#' \dontrun{
#' sim_gs_n(
#'   n_sim = 3,
#'   sample_size = 400,
#'   enroll_rate = enroll_rate,
#'   fail_rate = fail_rate,
#'   test = list(ia1 = ia1_test, ia2 = ia2_test, fa = fa_test),
#'   cut = list(ia1 = ia1_cut, ia2 = ia2_cut, fa = fa_cut)
#' )
#' }
#'
#' # Example 9: regular logrank test at all 3 analyses in parallel
#' plan("multisession", workers = 2)
#' sim_gs_n(
#'   n_sim = 3,
#'   sample_size = 400,
#'   enroll_rate = enroll_rate,
#'   fail_rate = fail_rate,
#'   test = wlr,
#'   cut = list(ia1 = ia1_cut, ia2 = ia2_cut, fa = fa_cut),
#'   weight = fh(rho = 0, gamma = 0)
#' )
#' plan("sequential")
#'
#' # Example 10: group sequential design with updated bounds
#' x <- gs_design_ahr(analysis_time = 1:3*12) |> to_integer()
#' sim_gs_n(
#'   n_sim = 3,
#'   sample_size = max(x$analysis$n),
#'   enroll_rate = x$enroll_rate,
#'   fail_rate = x$fail_rate,
#'   test = wlr,
#'   cut = list(ia1 = ia1_cut, ia2 = ia2_cut, fa = fa_cut),
#'   weight = fh(rho = 0, gamma = 0),
#'   original_design = x
#' )
sim_gs_n <- function(
    n_sim = 1000,
    sample_size = 500,
    stratum = data.frame(stratum = "All", p = 1),
    enroll_rate = data.frame(duration = c(2, 2, 10), rate = c(3, 6, 9)),
    fail_rate = data.frame(
      stratum = "All",
      duration = c(3, 100),
      fail_rate = log(2) / c(9, 18),
      hr = c(.9, .6),
      dropout_rate = rep(.001, 2)
    ),
    block = rep(c("experimental", "control"), 2),
    test = wlr,
    cut = NULL,
    original_design = NULL,
    ia_alpha_spending = c("min_planned_actual", "actual"),
    fa_alpha_spending = c("full_alpha", "info_frac"),
    ...) {
  # Input checking
  # TODO
  ia_alpha_spending <- match.arg(ia_alpha_spending)
  fa_alpha_spending <- match.arg(fa_alpha_spending)

  # parallel computation message for backends ----
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

  # Simulate for `n_sim` times
  ans <- foreach::foreach(
    sim_id = seq_len(n_sim),
    test = replicate(n=n_sim, expr=test, simplify = FALSE),
    .errorhandling = "stop",
    .options.future = list(seed = TRUE)
  ) %dofuture% {
    # Generate data
    simu_data <- sim_pw_surv(
      n = sample_size,
      stratum = stratum,
      block = block,
      enroll_rate = enroll_rate,
      fail_rate = to_sim_pw_surv(fail_rate)$fail_rate,
      dropout_rate = to_sim_pw_surv(fail_rate)$dropout_rate
    )

    # Initialize the cut date of IA(s) and FA
    n_analysis <- length(cut)
    cut_date <- rep(-100, n_analysis)
    ans_1sim <- NULL
    event_tbl <- NULL
    # Organize tests for each cutting
    if (is.function(test)) {
      test_single <- test
      test <- vector(mode = "list", length = n_analysis)
      test[] <- list(test_single)
    }
    if (length(test) != length(cut)) {
      stop("If you want to run different tests at each cutting, the list of
           tests must be the same length as the list of cuttings")
    }

    for (i_analysis in seq_len(n_analysis)) {
      # Get cut date
      cut_date[i_analysis] <- cut[[i_analysis]](simu_data)

      # Cut the data
      simu_data_cut <- simu_data |> cut_data_by_date(cut_date[i_analysis])

      # Test
      ans_1sim_new <- test[[i_analysis]](simu_data_cut, ...)
      ans_1sim_new <- c(sim_id = sim_id, ans_1sim_new)
      ans_1sim_new <- append(
        x = ans_1sim_new,
        values = c(
          analysis = i_analysis,
          cut_date = cut_date[i_analysis],
          n = nrow(simu_data_cut),
          event = sum(simu_data_cut$event)
        ),
        after = 3
      )
      ans_1sim_new <- convert_list_to_df_w_list_cols(ans_1sim_new)

      # bind rows of simulation results for all IA(s) and FA in 1 simulation
      ans_1sim_list <- list(ans_1sim, ans_1sim_new)
      ans_1sim <- rbindlist(ans_1sim_list, use.names = TRUE, fill = TRUE)

      # Get event counts per piecewise interval
      if (!is.null(original_design)){
        pw_event <- sapply(cumsum(original_design$fail_rate$duration),
                           function(threshold) {sum(simu_data_cut$tte < threshold)})
        event_tbl_new <- data.frame(analysis = rep(i_analysis, length(pw_event)),
                                    event = diff(c(0, pw_event)))
        event_tbl <- rbind(event_tbl, event_tbl_new)
      }
    }

    # Add planned and updated bounds
    if (!is.null(original_design)){
      # Add planned bounds
      planed_upper_bound <- original_design$z[original_design$bound == "upper"]
      planed_lower_bound <- original_design$z[original_design$bound == "lower"]
      ans_1sim$planed_upper_bound <- planed_upper_bound
      ans_1sim$planed_lower_bound <- planed_lower_bound

      # Calculate ustime and lstime
      obs_event <- (event_tbl |> dplyr::group_by(analysis) |> dplyr::summarize(x = sum(event)))$x
      plan_event <- original_design$analysis$event
      if (ia_alpha_spending == "actual" && fa_alpha_spending == "info_frac"){
        ustime <- obs_event / plan_event[n_analysis]
      } else if (ia_alpha_spending == "actual" && fa_alpha_spending == "full_alpha"){
        ustime <- obs_event / last(plan_event)
        ustime[n_analysis] <- 1
      } else if (ia_alpha_spending == "min_planned_actual" && fa_alpha_spending == "info_frac") {
        ustime <- pmin(obs_event, plan_event) / plan_event[n_analysis]
      } else {
        ustime <- pmin(obs_event, plan_event) / plan_event[n_analysis]
        ustime[n_analysis] <- 1
      }
      lstime <- ustime

      # Add updated bounds
      updated_design <- gsDesign2::gs_update_ahr(x = original_design,
                                                 alpha = original_design$input$alpha,
                                                 ustime = ustime, lstime = lstime,
                                                 event_tbl = event_tbl)$bound
      updated_upper_bound <- updated_design$z[updated_design$bound == "upper"]
      updated_lower_bound <- updated_design$z[updated_design$bound == "lower"]
      ans_1sim$updated_upper_bound <- updated_upper_bound
      ans_1sim$updated_lower_bound <- updated_lower_bound
    }

    ans_1sim
  }

  ans <- rbindlist(ans)
  setDF(ans)

  test_method <- ans$method[ans$sim_id == 1]
  if (all(substr(test_method, 1, 3) == "WLR")) {
    class(ans) <- c("simtrial_gs_wlr", class(ans))
    attr(ans, "method") <- unique(ans$parameter[ans$sim_id == 1])
  }

  return(ans)
}

#' Create a cutting function
#'
#' Create a cutting function for use with [sim_gs_n()]
#'
#' @param ... Arguments passed to [get_analysis_date()]
#'
#' @return A function that accepts a data frame of simulated trial data and
#'   returns a cut date
#'
#' @export
#'
#' @seealso [get_analysis_date()], [sim_gs_n()]
#'
#' @examples
#' # Simulate trial data
#' trial_data <- sim_pw_surv()
#'
#' # Create a cutting function that applies the following 2 conditions:
#' # - At least 45 months have passed since the start of the study
#' # - At least 300 events have occurred
#' cutting <- create_cut(
#'   planned_calendar_time = 45,
#'   target_event_overall = 350
#' )
#'
#' # Cut the trial data
#' cutting(trial_data)
create_cut <- function(...) {
  # Force evaluation of input arguments (required for parallel computing)
  lapply(X = list(...), FUN = force)

  function(data) {
    get_analysis_date(data, ...)
  }
}

#' Create a cutting test function
#'
#' Create a cutting test function for use with [sim_gs_n()]
#'
#' @param test A test function such as [wlr()], [maxcombo()], or [rmst()]
#' @param ... Arguments passed to the cutting test function
#'
#' @return A function that accepts a data frame of simulated trial data and
#'   returns a test result
#'
#' @export
#'
#' @seealso [sim_gs_n()], [create_cut()]
#'
#' @examples
#' # Simulate trial data
#' trial_data <- sim_pw_surv()
#'
#' # Cut after 150 events
#' trial_data_cut <- cut_data_by_event(trial_data, 150)
#'
#' # Create a cutting test function that can be used by sim_gs_n()
#' regular_logrank_test <- create_test(wlr, weight = fh(rho = 0, gamma = 0))
#'
#' # Test the cutting
#' regular_logrank_test(trial_data_cut)
#'
#' # The results are the same as directly calling the function
#' stopifnot(all.equal(
#'   regular_logrank_test(trial_data_cut),
#'   wlr(trial_data_cut, weight = fh(rho = 0, gamma = 0))
#' ))
create_test <- function(test, ...) {
  stopifnot(is.function(test))
  function(data) {
    test(data, ...)
  }
}

#' Perform multiple tests on trial data cutting
#'
#' WARNING: This experimental function is a work-in-progress. The function
#' arguments and/or returned output format may change as we add additional
#' features.
#'
#' @param data Trial data cut by [cut_data_by_event()] or [cut_data_by_date()]
#' @param ... One or more test functions. Use [create_test()] to change
#'   the default arguments of each test function.
#'
#' @return A list of test results, one per test. If the test functions are named
#'   in the call to `multitest()`, the returned list uses the same names.
#'
#' @export
#'
#' @seealso [create_test()]
#'
#' @examples
#' trial_data <- sim_pw_surv(n = 200)
#' trial_data_cut <- cut_data_by_event(trial_data, 150)
#'
#' # create cutting test functions
#' wlr_partial <- create_test(wlr, weight = fh(rho = 0, gamma = 0))
#' rmst_partial <- create_test(rmst, tau = 20)
#' maxcombo_partial <- create_test(maxcombo, rho = c(0, 0), gamma = c(0, 0.5))
#'
#' multitest(
#'   data = trial_data_cut,
#'   wlr = wlr_partial,
#'   rmst = rmst_partial,
#'   maxcombo = maxcombo_partial
#' )
multitest <- function(data, ...) {
  tests <- list(...)
  output <- vector(mode = "list", length = length(tests))
  names(output) <- names(tests)
  for (i in seq_along(tests)) {
    output[[i]] <- tests[[i]](data)
  }
  return(output)
}

# Convert a list to a one row data frame using list columns
convert_list_to_df_w_list_cols <- function(x) {
  stopifnot(is.list(x), !is.data.frame(x))

  new_list <- vector(mode = "list", length = length(x))
  names(new_list) <- names(x)

  for (i in seq_along(x)) {
    if (length(x[[i]]) > 1) {
      new_list[[i]] <- I(list(x[[i]]))
    } else {
      new_list[i] <- x[i]
    }
  }

  # Convert the list to a data frame with one row
  df_w_list_cols <- do.call(data.frame, new_list)
  stopifnot(nrow(df_w_list_cols) == 1)

  return(df_w_list_cols)
}
