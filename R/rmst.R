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

#' RMST difference of 2 arms
#'
#' @param data A time-to-event dataset with a column `tte` indicating the
#'   survival time and a column of `event` indicating whether it is
#'   event or censor.
#' @param tau RMST analysis time.
#' @param var_label_tte Column name of the TTE variable.
#' @param var_label_event Column name of the event variable.
#' @param var_label_group Column name of the grouping variable.
#' @param formula (default: `NULL`) A formula that indicates the TTE, event, and
#'   group variables using the syntax `Surv(tte, event) ~ group)` (see Details
#'   below for more information). This is an alternative to specifying the
#'   variables as strings. If a formula is provided, the values passed to
#'   `var_label_tte`, `var_label_event`, and `var_label_group` are ignored.
#' @param reference A group label indicating the reference group.
#' @param alpha Type I error.
#'
#' @details
#' The argument `formula` is provided as a convenience to easily specify the
#' TTE, event, and grouping variables using the syntax `Surv(tte, event) ~
#' group)`. `Surv()` is from the \{survival\} package ([survival::Surv()]).
#' You can also explicitly
#' name the arguments passed to `Surv()`, for example the following is
#' equivalent `Surv(event = event, time = tte) ~ group)`. Note however that the
#' function `Surv()` is never actually executed. Similarly, any other functions
#' applied in the formula are also ignored, thus you shouldn't apply any
#' transformation functions such as `log()` since these will have no effect.
#'
#' @return The z statistics.
#'
#' @export
#'
#' @examples
#' data(ex1_delayed_effect)
#' rmst(
#'   data = ex1_delayed_effect,
#'   var_label_tte = "month",
#'   var_label_event = "evntd",
#'   var_label_group = "trt",
#'   tau = 10,
#'   reference = "0"
#' )
#'
#' # Formula interface
#' rmst(
#'   data = ex1_delayed_effect,
#'   formula = Surv(month, evntd) ~ trt,
#'   tau = 10,
#'   reference = "0"
#' )
rmst <- function(
    data,
    tau = 10,
    var_label_tte = "tte",
    var_label_event = "event",
    var_label_group = "treatment",
    formula = NULL,
    reference = "control",
    alpha = 0.05) {
  stopifnot(is.data.frame(data))

  if (!is.null(formula)) {
    stopifnot(inherits(formula, "formula"))
    variables <- colnames(stats::get_all_vars(formula = formula, data = data))
    if (length(variables) != 3) {
      stop("The formula interface requires exactly 3 variables specified")
    }
    var_label_tte <- variables[1]
    var_label_event <- variables[2]
    var_label_group <- variables[3]
    # Properly parse the formula
    variables_parsed <- parse_formula_rmst(formula)
    var_label_tte <- variables_parsed["var_label_tte"]
    var_label_event <- variables_parsed["var_label_event"]
    var_label_group <- variables_parsed["var_label_group"]
  }

  res <- rmst_two_arm(
    time_var = data[[var_label_tte]],
    event_var = data[[var_label_event]],
    group_var = data[[var_label_group]],
    trunc_time = tau,
    reference = reference,
    alpha = alpha
  )

  ans <- list()
  ans$method <- "RMST"
  ans$parameter <- tau
  ans$estimate <- res$rmst_diff$rmst_diff
  ans$se <- res$rmst_diff$std
  ans$z <- res$rmst_diff$rmst_diff / res$rmst_diff$std

  return(ans)
}

#' Calculate RMST difference
#'
#' @inheritParams rmst_single_arm
#' @param group_var A vector of treatment groups.
#' @param trunc_time A numeric vector of pre-defined cut-off time point(s).
#' @param reference Group name of reference group for RMST comparison.
#'   Default is the first group name by alphabetical order.
#'
#' @return
#' A list of 2 data frames of RMST calculations:
#' - `rmst_per_arm`: the calculation results per group.
#' - `rmst_diff`: the calculation results of RMST differences.
#'
#' @importFrom data.table setDF rbindlist
#' @keywords internal
#'
#' @examples
#' data(ex1_delayed_effect)
#' with(
#'   ex1_delayed_effect,
#'   simtrial:::rmst_two_arm(
#'     time_var = month,
#'     event_var = evntd,
#'     group_var = trt,
#'     trunc_time = 6,
#'     reference = "0",
#'     alpha = 0.05
#'   )
#' )
rmst_two_arm <- function(
    time_var,
    event_var,
    group_var,
    trunc_time,
    reference = sort(unique(group_var))[1],
    alpha = 0.05) {
  # Input type check
  check_args(time_var, type = c("integer", "numeric"))
  check_args(event_var, type = c("integer", "numeric"))
  check_args(trunc_time, type = c("integer", "numeric"))
  check_args(reference, type = c("character"))
  check_args(alpha, type = c("integer", "numeric"))

  # Input value check
  stopifnot(time_var >= 0)
  stopifnot(event_var %in% c(0, 1))
  stopifnot(trunc_time >= 0)
  stopifnot(0 <= alpha & alpha <= 1)

  # Check truncation time
  if (any(min(tapply(time_var, group_var, max)) < trunc_time)) {
    stop(paste0(
      "The truncation time must be shorter than the minimum of the largest observed time in each group: ",
      sprintf("%.3f", min(tapply(time_var, group_var, max)))
    ))
  }

  g_label <- sort(unique(as.character(group_var)))

  op_single_list <- lapply(
    X = g_label,
    FUN = one_rmst,
    # arguments passed to FUN (one_rmst)
    group_var = group_var,
    time_var = time_var,
    event_var = event_var,
    trunc_time = trunc_time,
    alpha = alpha
  )
  op_single <- rbindlist(op_single_list)
  setDF(op_single)

  op_diff_list <- lapply(
    X = setdiff(g_label, reference),
    FUN = diff_rmst,
    # arguments passed to FUN (diff_rmst)
    op_single = op_single,
    reference = reference,
    trunc_time = trunc_time,
    alpha = alpha
  )
  op_diff <- rbindlist(op_diff_list)
  setDF(op_diff)

  ans <- list(rmst_per_arm = op_single, rmst_diff = op_diff)

  return(ans)
}

# Calculate RMST for each group by rmst_single()
one_rmst <- function(x, group_var, time_var, event_var, trunc_time, alpha) {
  indx <- group_var == x
  rmst_single_arm(
    time_var = time_var[indx],
    event_var = event_var[indx],
    tau = trunc_time,
    group_label = x,
    alpha = alpha
  )
}

# Calculate RMST difference and corresponding confidence intervals between each
# group with reference group
diff_rmst <- function(x, op_single, reference, trunc_time, alpha = alpha) {
  df_rf <- op_single[op_single$group == reference, ]
  df2 <- op_single[op_single$group == x, ]
  cutoff_time <- trunc_time
  group <- paste(unique(df2$group), "-", unique(df_rf$group))
  rmst_diff <- df2$rmst - df_rf$rmst
  variance <- df2$variance + df_rf$variance
  std <- sqrt(variance)
  lcl <- rmst_diff - stats::qnorm(1 - alpha / 2) * std
  ucl <- rmst_diff + stats::qnorm(1 - alpha / 2) * std
  data.frame(cutoff_time, group, rmst_diff, variance, std, lcl, ucl, stringsAsFactors = FALSE)
}

#' Calculate RMST for a single cut-off time point
#'
#' @param time_var A numeric vector of follow up time.
#' @param event_var A numeric or integer vector of the status indicator;
#'   0=alive 1=event.
#' @param tau A value of pre-defined cut-off time point.
#' @param group_label A character of customized treatment group name.
#' @param alpha A numeric value of the significant level for RMST
#'   confidence interval. Default is 0.05.
#'
#' @return
#' A data frame of
#' - Cutoff time: same as `tau`;
#' - Group label: same as `group_label`;
#' - Estimated RMST;
#' - Variance, standard error, and CIs of the estimated RMST;
#' - Number of events.
#'
#' @importFrom survival survfit Surv
#'
#' @keywords internal
#'
#' @examples
#' data(ex1_delayed_effect)
#' data_single_arm <- ex1_delayed_effect[ex1_delayed_effect$trt == 1, ]
#' simtrial:::rmst_single_arm(
#'   time_var = data_single_arm$month,
#'   event_var = data_single_arm$evntd,
#'   tau = 10,
#'   group_label = "Treatment 1",
#'   alpha = 0.05
#' )
rmst_single_arm <- function(
    time_var,
    event_var,
    tau,
    group_label = "Single Group",
    alpha = 0.05) {
  # Input type check
  check_args(time_var, type = c("integer", "numeric"))
  check_args(event_var, type = c("integer", "numeric"))
  check_args(tau, type = c("integer", "numeric"), length = 1)
  check_args(group_label, type = c("character", "factor"), length = 1)
  check_args(alpha, type = c("integer", "numeric"))

  # Input value check
  stopifnot(time_var >= 0)
  stopifnot(event_var %in% c(0, 1))
  stopifnot(tau >= 0)
  stopifnot(0 <= alpha & alpha <= 1)

  # Fit a single Kaplan-Meier curve
  fit <- survival::survfit(survival::Surv(time_var, event_var) ~ 1)

  # Extract survival probability, number of event, number at risk,
  # and number of censored along with observed time from the fitted model
  # as a new data frame.
  df <- data.frame(
    time = fit$time,
    n_risk = fit$n.risk,
    n_event = fit$n.event,
    n_censor = fit$n.censor,
    surv = fit$surv,
    stringsAsFactors = FALSE
  )

  # Filter df by survival time less or equal to the pre-specified cut-off time point tau
  df_fit1 <- df[df$time <= tau, ]

  # Add initial value of (time, survival) = (0,1) for calculating time difference
  df_fit2 <- rbind(df_fit1, c(0, NA, NA, NA, 1))

  # Add cut-off time if no records observed on the pre-specified time point
  if (max(df_fit1$time) != tau) df_fit2 <- rbind(df_fit2, c(tau, NA, NA, NA, NA))

  # Sort the data frame by time
  df_fit2 <- df_fit2[order(df_fit2$time), ]
  n_event <- df_fit2$n_event
  n_risk <- df_fit2$n_risk

  # Calculate the time difference and set the last value as NA
  time_diff <- c(diff((sort(df_fit2$time))), NA)

  # Calculate the area under the curve per time interval
  area <- time_diff * df_fit2$surv

  # Calculate the inverse cumulated area under the curve per observed time point A_i
  big_a <- rev(c(0, cumsum(rev(area)[-1])))

  # Calculation of dev refers to di / Yi * (Yi - di)
  dev <- (n_event / (n_risk * (n_risk - n_event))) * (big_a^2)

  # Based on the calculation, create a data frame with below items:
  # cutoff_time is the input of pre-defined cut-off time point
  cutoff_time <- tau
  # group is the input group name
  group <- group_label
  # rmst is the estimated RMST
  rmst <- sum(area, na.rm = TRUE)
  # std is the standard error of the estimated RMST
  variance <- sum(dev, na.rm = TRUE) * sum(n_event, na.rm = TRUE) / (sum(n_event, na.rm = TRUE) - 1)
  std <- sqrt(variance)
  # lcl and ucl are lower/upper control limit of CIs for RMST
  lcl <- rmst - stats::qnorm(1 - alpha / 2) * std
  ucl <- rmst + stats::qnorm(1 - alpha / 2) * std
  event <- sum(n_event, na.rm = TRUE)

  ans <- data.frame(
    cutoff_time, group, rmst, variance, std, lcl, ucl, event,
    stringsAsFactors = FALSE
  )

  return(ans)
}

#' Parse the formula argument of rmst()
#'
#' The canonical syntax is `Surv(tte, event) ~ group)`
#'
#' @inheritParams rmst
#'
#' @examples
#' parse_formula_rmst(formula = Surv(tte, event) ~ group)
#'
#' parse_formula_rmst(formula = Surv(event = event, tte) ~ group)
#'
#' @noRd
parse_formula_rmst <- function(formula) {
  formula_terms <- stats::terms(formula)

  var_label_group <- attr(formula_terms, "term.labels")[1]
  if (length(var_label_group) != 1) {
    stop(
      "Unable to identify a single group variable. Received: '",
      paste(var_label_group, collapse = " "), "'"
    )
  }

  surv_call <- attr(formula_terms, "variables")[[2]]
  if (is.symbol(surv_call)) {
    stop("Must use canonical formula syntax with Surv()")
  }
  surv_call_match <- match.call(definition = function(time, event) {}, surv_call)

  var_label_tte <- as.character(surv_call_match$time)
  if (length(var_label_tte) != 1) {
    stop(
      "Unable to identify a single tte variable. Received: '",
      paste(var_label_tte, collapse = " "), "'"
    )
  }
  var_label_event <- as.character(surv_call_match$event)
  if (length(var_label_event) != 1) {
    stop(
      "Unable to identify a single event variable. Received: '",
      paste(var_label_event, collapse = " "), "'"
    )
  }

  result <- c(
    "var_label_tte" = var_label_tte,
    "var_label_event" = var_label_event,
    "var_label_group" = var_label_group
  )

  return(result)
}
