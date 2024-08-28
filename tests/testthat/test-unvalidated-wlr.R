test_that("wlr() accepts tte_data and counting_process objects as input", {
  # cut_data_by_event()
  x <- sim_pw_surv(n = 300) |> cut_data_by_event(100)
  expect_s3_class(x, "tte_data")
  results_tte_data <- x |> wlr(weight = fh(0, 0.5))

  x <- x |> counting_process(arm = "experimental")
  expect_s3_class(x, "counting_process")
  results_counting_process <- x |> wlr(weight = fh(0, 0.5))

  expect_equal(results_tte_data, results_counting_process)

  # cut_data_by_date()
  x <- sim_pw_surv(n = 300) |> cut_data_by_date(cut_date = 300)
  expect_s3_class(x, "tte_data")
  results_tte_data <- x |> wlr(weight = fh(0, 0.5))

  x <- x |> counting_process(arm = "experimental")
  expect_s3_class(x, "counting_process")
  results_counting_process <- x |> wlr(weight = fh(0, 0.5))

  expect_equal(results_tte_data, results_counting_process)
})

test_that("wlr() rejects input object without proper class", {
  x <- mtcars
  expect_error(wlr(x), "no applicable method")
})
