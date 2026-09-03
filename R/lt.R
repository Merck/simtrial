#  Copyright (c) 2026 Merck & Co., Inc., Rahway, NJ, USA and its affiliates.
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

# Re-export the lt() generic so users can call lt() on a simulation summary after
# only loading simtrial (without also attaching lt or qualifying with lt::).
# This also makes S3 dispatch robust regardless of package load order.
# See https://github.com/yihui/lt/issues/4.

#' @importFrom lt lt
#' @export
lt::lt

#' Create an lt table from a simulation summary
#'
#' S3 method for [lt()] that converts a group sequential simulation summary
#' (a `simtrial_gs_wlr` object returned by [summary()]) into a formatted lt
#' table. This is the lightweight replacement for the deprecated [as_gt()].
#'
#' @param data A summary object returned by [summary()].
#' @param title Title of the lt table.
#' @param subtitle Subtitle of the lt table.
#' @param ... Additional arguments (not used).
#'
#' @return An `lt_tbl` object summarizing the simulation results.
#'
#' @name lt-methods
#'
#' @seealso [as_gt()]
#'
#' @exportS3Method lt::lt
#'
#' @examples
#'
#' # Parameters for enrollment
#' enroll_rampup_duration <- 4 # Duration for enrollment ramp up
#' enroll_duration <- 16 # Total enrollment duration
#' enroll_rate <- gsDesign2::define_enroll_rate(
#'   duration = c(
#'     enroll_rampup_duration, enroll_duration - enroll_rampup_duration),
#'  rate = c(10, 30))
#'
#' # Parameters for treatment effect
#' delay_effect_duration <- 3 # Delay treatment effect in months
#' median_ctrl <- 9 # Survival median of the control arm
#' median_exp <- c(9, 14) # Survival median of the experimental arm
#' dropout_rate <- 0.001
#' fail_rate <- gsDesign2::define_fail_rate(
#'   duration = c(delay_effect_duration, 100),
#'   fail_rate = log(2) / median_ctrl,
#'   hr = median_ctrl / median_exp,
#'   dropout_rate = dropout_rate)
#'
#' # Other related parameters
#' alpha <- 0.025 # Type I error
#' beta <- 0.1 # Type II error
#' ratio <- 1 # Randomization ratio (experimental:control)
#'
#' # Build a one-sided group sequential design
#' design <- gsDesign2::gs_design_ahr(
#'   enroll_rate = enroll_rate, fail_rate = fail_rate,
#'   ratio = ratio, alpha = alpha, beta = beta,
#'   analysis_time = c(12, 24, 36),
#'   upper = gsDesign2::gs_spending_bound,
#'   upar = list(sf = gsDesign::sfLDOF, total_spend = alpha),
#'   lower = gsDesign2::gs_b,
#'   lpar = rep(-Inf, 3))
#'
#' # Define cuttings of 2 IAs and 1 FA
#' ia1_cut <- create_cut(target_event_overall = ceiling(design$analysis$event[1]))
#' ia2_cut <- create_cut(target_event_overall = ceiling(design$analysis$event[2]))
#' fa_cut <- create_cut(target_event_overall = ceiling(design$analysis$event[3]))
#'
#' # Run simulations
#' simulation <- sim_gs_n(
#'   n_sim = 3,
#'   sample_size = ceiling(design$analysis$n[3]),
#'   enroll_rate = design$enroll_rate,
#'   fail_rate = design$fail_rate,
#'   test = wlr,
#'   cut = list(ia1 = ia1_cut, ia2 = ia2_cut, fa = fa_cut),
#'   weight = fh(rho = 0, gamma = 0.5))
#'
#' # Summarize simulations
#' simulation |>
#'  summary(bound = gsDesign::gsDesign(k = 3, test.type = 1, sfu = gsDesign::sfLDOF)$upper$bound) |>
#'  lt()
#'
#' # Summarize simulations and compare with the planned design
#' simulation |>
#'   summary(design = design) |>
#'   lt()
lt.simtrial_gs_wlr <- function(data,
                               title = "Summary of simulation results by WLR tests",
                               subtitle = NULL, ...){
  x <- data

  # The raw output of sim_gs_n() also carries the "simtrial_gs_wlr" class but is
  # not a summary (it lacks the attributes added by summary()). In that case fall
  # back to a plain lt table, mirroring how bare gt() used to render it.
  if (is.null(attributes(x)$compare_with_design)) {
    return(lt::lt(as.data.frame(x), ...))
  }

  # get the default subtitle
  if (is.null(subtitle)) {
    subtitle <- paste0("Weighted by ", attributes(x)$method)
  }

  # if it is not compared with the design
  if (attributes(x)$compare_with_design == "no") {
    as.data.frame(x) |>
      lt::lt() |>
      lt::lt_label(sim_time = "Time", sim_n = "N", sim_event = "Event", sim_upper_prob = "Crossing probability") |>
      lt::lt_move(columns = c("sim_time", "sim_n", "sim_event"), after = "analysis") |>
      lt::lt_header(title = title, subtitle = subtitle)
  } else {
    # get the design type, either one-sided or two-sided
    design_type <- attributes(x)$design_type

    # lt has no tidyselect, so enumerate the columns of each spanner explicitly.
    # The columns must be listed in the same paired order (asymptotic before
    # simulated) that lt_move() lays them out below, because lt matches a
    # spanner to the visual position of its first column and then spans the
    # next length(columns) columns; listing them in any other order would
    # misalign the spanners (and silently drop later ones).
    time_cols <- c("asy_time", "sim_time")
    n_cols <- c("asy_n", "sim_n")
    event_cols <- c("asy_event", "sim_event")
    upper_cols <- c("asy_upper_prob", "sim_upper_prob")
    lower_cols <- c("asy_lower_prob", "sim_lower_prob")

    # build an lt table as return, moving the paired asymptotic/simulated columns
    # right after `analysis` so each spanner covers a contiguous block
    ans <- as.data.frame(x) |>
      lt::lt() |>
      lt::lt_move(
        columns = c(time_cols, n_cols, event_cols),
        after = "analysis")

    # for a two-sided design, keep the efficacy (upper) and futility (lower)
    # probability columns contiguous within their own spanners
    if (design_type == "two-sided") {
      ans <- ans |>
        lt::lt_move(
          columns = c(upper_cols, lower_cols),
          after = "sim_event")
    }

    ans <- ans |>
      lt::lt_spanner(label = "Time", columns = time_cols) |>
      lt::lt_spanner(label = "Events", columns = event_cols) |>
      lt::lt_spanner(label = "N", columns = n_cols) |>
      lt::lt_spanner(
        label = "Probability of crossing efficacy bounds under H1",
        columns = upper_cols)

    if (design_type == "two-sided") {
      ans <- ans |> lt::lt_spanner(
        label = "Probability of crossing futility bounds under H1",
        columns = lower_cols)
    }

    # label the asymptotic/simulated/analysis columns, mirroring the
    # starts_with()/matches() rules used by gt::cols_label() in as_gt()
    labels <- ifelse(
      startsWith(names(x), "asy"), "Asymptotic",
      ifelse(startsWith(names(x), "sim"), "Simulated", "Analysis"))
    labels <- stats::setNames(as.list(labels), names(x))

    do.call(lt::lt_label, c(list(ans), labels)) |>
      lt::lt_header(title = title, subtitle = subtitle)
  }
}
