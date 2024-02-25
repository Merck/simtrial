# Helper functions used by test-independent_test_counting_process.R

surv_to_count <- function(time, status, trt, strats) {
  db <- data.frame(time, status, trt, strats)

  # KM estimator by stratum
  tidy_survfit <- function(...) {
    .survfit <- summary(survival::survfit(...))
    n <- length(.survfit$time)
    data.frame(
      time = .survfit$time,
      n.risk = .survfit$n.risk,
      n.event = .survfit$n.event,
      surv = c(1, .survfit$surv[-n])
    ) # ensure left continuous
  }
  km <- db |>
    dplyr::group_by(strats) |>
    dplyr::do(tidy_survfit(Surv(time, status) ~ 1, data = .))

  # KM estimator by stratum and treatment Group Predicted at Specified Time
  pred_survfit <- function(pred_time, ...) {
    .survfit <- survival::survfit(...)

    # At risk subjects at pred_time
    n.risk <- stepfun(.survfit$time, c(.survfit$n.risk, 0), right = TRUE)(pred_time)
    .x1 <- data.frame(time = pred_time, n.risk)

    # Number of Event
    .x2 <- data.frame(time = .survfit$time, n.event = .survfit$n.event) |> subset(n.event > 0)

    merge(.x1, .x2, all = TRUE) |> dplyr::mutate(n.event = dplyr::if_else(is.na(n.event), 0, n.event))
  }

  km_by_trt <- db |>
    dplyr::group_by(strats, trt) |>
    dplyr::do(pred_time = pred_survfit(km[km$strats == .$strats[1], ]$time,
      Surv(time, status) ~ 1,
      data = .
    )) |>
    tidyr::unnest(cols = pred_time) |>
    dplyr::rename(tn.risk = n.risk, tn.event = n.event)

  # Log Rank Expectation Difference and Variance
  res <- merge(km, km_by_trt, all = TRUE) |>
    dplyr::arrange(trt, strats, time) |>
    dplyr::mutate(
      OminusE = tn.event - tn.risk / n.risk * n.event,
      Var = (n.risk - tn.risk) * tn.risk * n.event * (n.risk - n.event) / n.risk^2 / (n.risk - 1)
    )
}
