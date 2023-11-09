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

#' Zero early weight for weighted logrank tests
#'
#' @param x A [counting_process()]-class `tibble` with a counting process dataset.
#' @param early_period The initial delay period where weights increase;
#'   after this, weights are constant at the final weight in the delay period.
#' @param fail_rate a data frame record the failure rate
#'
#' @return A data frame. The column `weight` contains the weights for the
#'   early zero weighted logrank test for the data in `x`.
#'
#' @references
#' Xu, Z., Zhen, B., Park, Y., & Zhu, B. (2017).
#' "Designing therapeutic cancer vaccine trials with delayed treatment effect."
#' _Statistics in medicine_, 36(4), 592-605.
#'
#' @export
#'
#' @examples
#' library(dplyr)
#' library(gsDesign2)
#'
#' # Example 1: unstratified
#' sim_pw_surv(n = 200) %>%
#'   cut_data_by_event(125) %>%
#'   counting_process(arm = "experimental") %>%
#'   early_zero_weight(early_period = 2) %>%
#'   filter(row_number() %in% seq(5, 200, 40))
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
#' enroll_rate$rate <- enroll_rate$rate * n / sum(enroll_rate$duration * enroll_rate$rate)
#'
#' # failure rate
#' med_pos <- 10         # median of the biomarker positive population
#' med_neg <- 8          # median of the biomarker negative population
#' hr_pos <- c(1, 0.7)   # hazard ratio of the biomarker positive population
#' hr_neg <- c(1, 0.8)   # hazard ratio of the biomarker negative population
#' fail_rate <- define_fail_rate(
#'   stratum = rep(stratum, each = 2),
#'   duration = c(3, 1000, 4, 1000),
#'   fail_rate = c(log(2) / c(med_pos, med_pos, med_neg, med_neg)),
#'   hr = c(hr_pos, hr_neg),
#'   dropout_rate = 0.01)
#'
#' # simulate data
#' temp <- simfix2simpwsurv(fail_rate)        # transfer the failure rate
#' set.seed(2023)
#'
#' sim_pw_surv(
#'   n = n,                                                            # sample size
#'   stratum = tibble(stratum = stratum, p = prevelance_ratio),        # stratified design with prevalence ratio of 6:4
#'   block =  c("control", "control", "experimental", "experimental"), # randomization ratio
#'   enroll_rate = enroll_rate,                                        # enrollment rate
#'   fail_rate = temp$fail_rate,                                       # failure rate
#'   dropout_rate = temp$dropout_rate                                  # dropout rate
#'   ) %>%
#'   cut_data_by_event(125) %>%
#'   counting_process(arm = "experimental") %>%
#'   early_zero_weight(early_period = 2, fail_rate = fail_rate) %>%
#'   filter(row_number() %in% seq(5, 200, 40))
early_zero_weight <- function(x, early_period = 4, fail_rate = NULL) {

  n_stratum <- length(unique((x$stratum)))

  # if it is unstratified design
  if(n_stratum == 1){
    ans <- x %>%
      mutate(weight = case_when(tte < early_period ~ 0,
                                tte >= early_period ~ 1))
  } else {
    if(is.null(fail_rate)){
      stop("For stratified design to enable `early_zero_weight`, users must input the `fail_rate`.")
    }
    if(!all((fail_rate %>% group_by(stratum) %>% summarise(x = n() == 2))$x)){
      stop("`early_zero_weight` only allow delayed treatment effect, i.e., 2 piece failure rate with HR = 1 at the first period.")
    }

    late_hr <- fail_rate %>%
      filter(hr != 1) %>%
      select(stratum, hr)
    delay_change_time <- fail_rate %>%
      filter(hr == 1) %>%
      select(stratum, duration)

    ans <- x %>%
      left_join(late_hr) %>%
      left_join(delay_change_time) %>%
      mutate(weight = case_when(tte < duration ~ 0,
                                tte >= duration ~ hr))

  }

  return(ans)
}
