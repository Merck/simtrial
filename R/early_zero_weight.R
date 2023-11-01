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

#' Magirr and Burman modestly weighted logrank tests
#'
#' Magirr and Burman (2019) proposed a weighted logrank test to have better
#' power than the logrank test when the treatment effect is delayed,
#' but to still maintain good power under a proportional hazards assumption.
#' In Magirr (2021), (the equivalent of) a maximum weight was proposed
#' as opposed to a fixed time duration over which weights would increase.
#' The weights for some early interval specified by the user are the inverse
#' of the combined treatment group empirical survival distribution; see details.
#' After this initial period, weights are constant at the maximum of the
#' previous weights. Another advantage of the test is that under strong
#' null hypothesis that the underlying survival in the control group is
#' greater than or equal to underlying survival in the experimental group,
#' Type I error is controlled as the specified level.
#'
#' Computes Magirr-Burman weights and adds them to a dataset created by
#' [counting_process()].
#' These weights can then be used to compute a z-statistic for the
#' modestly weighted logrank test proposed.
#'
#' @param x A [counting_process()]-class `tibble` with a counting process dataset.
#' @param delay The initial delay period where weights increase;
#'   after this, weights are constant at the final weight in the delay period.
#' @param w_max Maximum weight to be returned.
#'   Set `delay = Inf`, `w_max = 2` to be consistent with recommendation of
#'   Magirr (2021).
#'
#' @return A data frame. The column `mb_weight` contains the weights for the
#'   Magirr-Burman weighted logrank test for the data in `x`.
#'
#' @details
#' We define \eqn{t^*} to be the input variable `delay`.
#' This specifies an initial period during which weights increase.
#' We also set a maximum weight \eqn{w_{\max}}.
#' To define specific weights, we let \eqn{S(t)} denote the Kaplan-Meier
#' survival estimate at time \eqn{t} for the combined data
#' (control plus experimental treatment groups).
#' The weight at time \eqn{t} is then defined as
#' \deqn{w(t)=\min(w_{\max}, S(\min(t, t^*))^{-1}).}
#'
#' @references
#' Magirr, Dominic, and Carl‐Fredrik Burman. 2019.
#' "Modestly weighted logrank tests."
#' _Statistics in Medicine_ 38 (20): 3782--3790.
#'
#' Magirr, Dominic. 2021.
#' "Non‐proportional hazards in immuno‐oncology: Is an old perspective needed?"
#' _Pharmaceutical Statistics_ 20 (3): 512--527.
#'
#' @export
#'
#' @examples
#' library(dplyr)
#' library(gsDesign2)
#'
#' # Example 1: unstratified
#'
#' # Example 2: stratified
#' n <- 500
#' # two strata
#' stratum <- c("Biomarker-positive", "Biomarker-negative")
#' prevelance_ratio <- c(0.6, 0.4)
#'
#' # enrollment rate
#' enroll_rate <- define_enroll_rate(
#'   stratum = rep(stratum, each = 2),
#'   duration = c(2, 10, 2, 10),
#'   rate =  c(c(1, 4) * prevelance_ratio[1], c(1, 4) * prevelance_ratio[2]))
#'enroll_rate$rate <- enroll_rate$rate * n / sum(enroll_rate$duration * enroll_rate$rate)
#'
#' # failure rate
#' med_pos <- 10         # median of the biomarker positive population
#' med_neg <- 8          # median of the biomarker negative population
#' hr_pos <- c(1, 0.7)   # hazard ratio of the biomarker positive population
#' hr_neg <- c(1, 0.8)   # hazard ratio of the biomarker negative population
#' fail_rate <- define_fail_rate(
#'   stratum = rep(stratum, each = 2),
#'   duration = 1000,
#'   fail_rate = c(log(2) / c(med_pos, med_pos, med_neg, med_neg)),
#'   hr = c(hr_pos, hr_neg),
#'   dropout_rate = 0.01)
#'
#' # simulate data
#' temp <- simfix2simpwsurv(fail_rate)        # transfer the failure rate
#' set.seed(2023)
#'
#' x <- sim_pw_surv(
#'   n = n,                                                            # sample size
#'   stratum = tibble(stratum = stratum, p = prevelance_ratio),        # stratified design with prevalence ratio of 6:4
#'   block =  c("control", "control", "experimental", "experimental"), # randomization ratio
#'   enroll_rate = enroll_rate,                                        # enrollment rate
#'   fail_rate = temp$fail_rate,                                       # failure rate
#'   dropout_rate = temp$dropout_rate                                  # dropout rate
#'   ) %>%
#'   cut_data_by_event(125) %>%
#'   counting_process(arm = "experimental")
early_zero_weight <- function(x, delay = 4, w_max = Inf) {

  # Compute max weight by stratum
  x2 <- x %>% group_by(stratum)
  # Make sure you don't lose any stratum!
  tbl_all_stratum <- x2 %>% summarize()

  x2 %>%
    # Look only up to delay time
    filter(tte <= 6) %>%
    # Weight before delay specified as 1/S
    summarize(max_weight = max(1 / s)) %>%
    # Get back stratum with no records before delay ends
    right_join(tbl_all_stratum, by = "stratum") %>%
    # `max_weight` is 1 when there are no records before delay ends
    mutate(max_weight = replace_na(max_weight, 1)) %>%
    # Cut off weights at w_max
    mutate(max_weight = pmin(Inf, max_weight)) %>%
    # Now merge max_weight back to stratified dataset
    full_join(x2, by = "stratum") %>%
    # Weight is min of max_weight and 1/S which will increase up to delay
    mutate(mb_weight = pmin(max_weight, 1 / s)) %>%
    select(-max_weight)


}
