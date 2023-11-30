# Prepare `mb_delayed_effect` dataset
library(simtrial)

# Load existing object for comparison
load("data/mb_delayed_effect.rda")
existing <- mb_delayed_effect

set.seed(6671)

ds <- sim_pw_surv(
  n = 200,
  block = c(rep("control", 2), rep("experimental", 2)),
  enroll_rate = data.frame(rate = 200 / 12, duration = 12),
  fail_rate = data.frame(
    stratum = c("All", "All", "All"),
    period = c(1, 1, 2),
    treatment = c("control", "experimental", "experimental"),
    duration = c(42, 6, 36),
    rate = c(log(2) / 15, log(2) / 15, log(2) / 15 * 0.6)
  ),
  dropout_rate = data.frame(
    stratum = c("All", "All"),
    period = c(1, 1),
    treatment = c("control", "experimental"),
    duration = c(42, 42),
    rate = c(0, 0)
  )
)

# Cut data at 24 months after final enrollment
mb_delayed_effect <- cut_data_by_date(ds, max(ds$enroll_time) + 24)

if (!all.equal(existing, mb_delayed_effect)) {
  warning(
    "The updated mb_delayed_effect differs from the existing object",
    .call = FALSE,
    immediate. = TRUE
  )
}

usethis::use_data(mb_delayed_effect, overwrite = TRUE)
