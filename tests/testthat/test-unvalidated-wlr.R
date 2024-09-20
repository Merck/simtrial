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

test_that("wlr() uses argument ratio", {
  x <- data.frame(treatment = ifelse(ex1_delayed_effect$trt == 1, "experimental", "control"),
                 stratum = rep("All", nrow(ex1_delayed_effect)),
                 tte = ex1_delayed_effect$month,
                 event = ex1_delayed_effect$evntd)
  class(x) <- c("tte_data", class(x))
  wlr_w_ratio <- x |> wlr(weight = fh(rho = 0, gamma = 0.5), ratio = 2)
  wlr_wo_ratio <- x |> wlr(weight = fh(rho = 0, gamma = 0.5))
  expect_false(isTRUE(all.equal(wlr_w_ratio, wlr_wo_ratio)))
})
