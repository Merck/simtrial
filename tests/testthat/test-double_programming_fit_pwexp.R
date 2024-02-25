test_that("fit_pwexp works when the input vector `intervals` contains a final infinite time interval", {
  Srv <- Surv(time = ex2_delayed_effect$month, event = ex2_delayed_effect$evntd)
  expect_equal(fit_pwexp(intervals = c(3, 6, Inf)), fit_pwexp(intervals = c(3, 6, Inf)))
})

test_that("fit_pwexp works when at least one event occurred after sum(intervals)", {
  Srv <- Surv(time = ex2_delayed_effect$month, event = ex2_delayed_effect$evntd)
  expect_equal(fit_pwexp(intervals = c(3, 6, 6)), fit_pwexp(intervals = c(3, 6, 6)))
})

test_that("fit_pwexp works when sum(intervals) covers all events in the observed data", {
  Srv <- Surv(time = ex2_delayed_effect$month, event = ex2_delayed_effect$evntd)
  expect_equal(fit_pwexp(intervals = c(3, 6, 50)), fit_pwexp(intervals = c(3, 6, 50)))
})

test_that("fit_pwexp works when no events observed in some pieces of time in the input vector `intervals`", {
  Srv <- Surv(time = ex2_delayed_effect$month, event = ex2_delayed_effect$evntd)
  max <- max(Srv[, "time"])
  expect_equal(fit_pwexp(intervals = c(max, max + 3, max + 3)), fit_pwexp(intervals = c(max, max + 3, max + 3)))
})
