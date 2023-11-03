## code to prepare `MBdelayed` dataset goes here
library(simtrial)
library(tibble)
set.seed(6671)
ds <- sim_pw_surv(
  n = 200,
  block = c(rep("control", 2), rep("experimental", 2)),
  enroll_rate = tibble(rate = 200 / 12, duration = 12),
  fail_rate = tribble(
    ~stratum, ~period, ~treatment, ~duration, ~rate,
    "All", 1, "control", 42, log(2) / 15,
    "All", 1, "experimental", 6, log(2) / 15,
    "All", 2, "experimental", 36, log(2) / 15 * 0.6
  ),
  dropout_rate = tribble(
    ~stratum, ~period, ~treatment, ~duration, ~rate,
    "All", 1, "control", 42, 0,
    "All", 1, "experimental", 42, 0
  )
)
# cut data at 24 months after final enrollment
MBdelayed <- ds %>% cut_data_by_date(max(ds$enroll_time) + 24)

usethis::use_data(MBdelayed, overwrite = TRUE)
