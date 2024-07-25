test_that("Verify the info output from WLR-early-zero test", {
  data <- sim_pw_surv(n = 200) |> cut_data_by_event(100)

  observed <- data |> wlr(weight = early_zero(early_period = 6))

  expected_info <- data |>
    filter(tte >= 6) |>
    group_by(treatment) |>
    summarize(event = sum(event)) |>
    ungroup() |>
    summarize(info = 1 / sum(1 / event))

  expected_info0 <- data |>
    filter(tte >= 6) |>
    summarize(event = sum(event)) |>
    summarize(info = event * 0.5 * 0.5)

  expect_equal(observed$info0, expected_info0)
  expect_equal(observed$info, expected_info)
})

