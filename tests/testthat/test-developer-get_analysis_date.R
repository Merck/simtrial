test_that("Cut by targeted event per stratum", {
  sim_data <- sim_pw_surv(n = 360,
                          stratum = data.frame(stratum = c("Positive", "Negative"), p = c(1/3, 2/3)),
                          block = c(rep("control", 2), rep("experimental", 2)),
                          enroll_rate = data.frame(stratum = c("Positive", "Negative"),
                                                   duration = c(12, 12),
                                                   rate = c(10, 20)),
                          fail_rate = data.frame(stratum = rep(c("Positive", "Negative"), each = 2),
                                                 period = c(1, 1, 1, 1),
                                                 treatment = c("control", "experimental", "control", "experimental"),
                                                 duration = c(Inf, Inf, Inf, Inf),
                                                 rate = log(2) / c(9, 12, 6, 10)),
                          dropout_rate = data.frame(stratum = rep(c("Positive", "Negative"), each = 2),
                                                    period = c(1, 1, 1, 1),
                                                    treatment = c("control", "experimental", "control", "experimental"),
                                                    duration = c(Inf, Inf, Inf, Inf),
                                                    rate = c(0.001, 0.001, 0.001, 0.001)
                          ))

  # test 1
  cut_date <- sim_data |> get_analysis_date(target_event_per_stratum = c("Positive" = 50, "Negative" = NA))
  cut_data <- sim_data |> cut_data_by_date(cut_date)
  expect_equal(50,
               sum(cut_data[which(cut_data$stratum == "Positive"), ]$event)
               )

  # test 2
  cut_date <- sim_data |> get_analysis_date(target_event_per_stratum = c("Positive" = NA, "Negative" = 50))
  cut_data <- sim_data |> cut_data_by_date(cut_date)
  expect_equal(50,
               sum(cut_data[which(cut_data$stratum == "Negative"), ]$event)
               )

  # test 3
  cut_date <- sim_data |> get_analysis_date(target_event_per_stratum = c("Positive" = 50, "Negative" = 70))
  cut_data <- sim_data |> cut_data_by_date(cut_date)
  expect_true(sum(cut_data[which(cut_data$stratum == "Positive"), ]$event) >= 50)
  expect_true(sum(cut_data[which(cut_data$stratum == "Negative"), ]$event) >= 70)
})

test_that("Cut by targeted sample size per stratum", {
  sim_data <- sim_pw_surv(n = 360,
                          stratum = data.frame(stratum = c("Positive", "Negative"), p = c(1/3, 2/3)),
                          block = c(rep("control", 2), rep("experimental", 2)),
                          enroll_rate = data.frame(stratum = c("Positive", "Negative"),
                                                   duration = c(12, 12),
                                                   rate = c(10, 20)),
                          fail_rate = data.frame(stratum = rep(c("Positive", "Negative"), each = 2),
                                                 period = c(1, 1, 1, 1),
                                                 treatment = c("control", "experimental", "control", "experimental"),
                                                 duration = c(Inf, Inf, Inf, Inf),
                                                 rate = log(2) / c(9, 12, 6, 10)),
                          dropout_rate = data.frame(stratum = rep(c("Positive", "Negative"), each = 2),
                                                    period = c(1, 1, 1, 1),
                                                    treatment = c("control", "experimental", "control", "experimental"),
                                                    duration = c(Inf, Inf, Inf, Inf),
                                                    rate = c(0.001, 0.001, 0.001, 0.001)
                          ))

  # test 1
  cut_date <- sim_data |> get_analysis_date(min_n_per_stratum = c("Positive" = 50, "Negative" = 100),
                                            min_followup = 6)
  cut_data <- sim_data |> cut_data_by_date(cut_date)
  expect_true(cut_data |> dplyr::filter(stratum == "Positive") |> nrow() >= 50)
  expect_true(cut_data |> dplyr::filter(stratum == "Negative") |> nrow() >= 100)

  # test 2
  cut_date <- sim_data |> get_analysis_date(min_n_per_stratum = c("Positive" = 50, "Negative" = NA),
                                            min_followup = 6)
  cut_data <- sim_data |> cut_data_by_date(cut_date)
  expect_true(cut_data |> dplyr::filter(stratum == "Positive") |> nrow() >= 50)

  # test 3
  cut_date <- sim_data |> get_analysis_date(target_event_per_stratum = c("Positive" = NA, "Negative" = 100))
  cut_data <- sim_data |> cut_data_by_date(cut_date)
  expect_true(cut_data |> dplyr::filter(stratum == "Negative") |> nrow() >= 100)
})

