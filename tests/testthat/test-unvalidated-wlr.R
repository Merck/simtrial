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

test_that("wlr() rejects input object without proper columns", {
  x <- mtcars
  expect_error(
    wlr(x),
    'Input must have the columns "tte", "event", and "treatment"'
  )
})

test_that("wlr() accepts unclassed input object with proper columns", {
  # Users should be able to pass unclassed custom objects
  x <- sim_pw_surv(n = 300) |> cut_data_by_event(100)
  expected <- wlr(x, weight = fh(0, 0.5))
  class(x) <- "data.frame"
  observed <- wlr(x, weight = fh(0, 0.5))
  expect_equal(observed, expected)
})

test_that("wlr() uses argument ratio", {
  x <- data.frame(
    treatment = ifelse(ex1_delayed_effect$trt == 1, "experimental", "control"),
    stratum = rep("All", nrow(ex1_delayed_effect)),
    tte = ex1_delayed_effect$month,
    event = ex1_delayed_effect$evntd
  )
  wlr_w_ratio <- x |> wlr(weight = fh(rho = 0, gamma = 0.5), ratio = 2)
  wlr_wo_ratio <- x |> wlr(weight = fh(rho = 0, gamma = 0.5))
  expect_false(isTRUE(all.equal(wlr_w_ratio, wlr_wo_ratio)))
})

test_that("cut_data_by_date() and cut_data_by_event() return the same classes", {
  x <- sim_pw_surv(n = 300)
  data_by_event <- cut_data_by_event(x, 100)
  data_by_date <- cut_data_by_date(x, cut_date = 300)

  expect_identical(class(data_by_event), class(data_by_date))
})

test_that("wlr() formula argument can rename columns", {
  x <- sim_pw_surv(n = 300) |> cut_data_by_event(100)
  expected <- wlr(x, weight = fh(0, 0.5))

  # Rearrange and rename columns to simulate custom user data. Also remove class
  # "tte_data"
  y <- x[, rev(colnames(x))]
  colnames(y) <- toupper(colnames(y))
  class(y) <- "data.frame"

  observed <- wlr(y, weight = fh(0, 0.5), formula =  Surv(TTE, EVENT) ~ TREATMENT)

  # Sometimes info0 is off by ~2e-5
  expect_equal(observed, expected, tolerance = 1e-4)
})

test_that("wlr() accepts formula for unstratified design", {
  data("ex1_delayed_effect", package = "simtrial", envir = environment())
  ex1_delayed_effect$trt <- ifelse(
    ex1_delayed_effect$trt == 1,
    "experimental",
    "control"
  )

  # Use ex1_delayed_effect directly via the argument `formula`
  observed <- wlr(
    data = ex1_delayed_effect,
    formula = Surv(month, evntd) ~ trt,
    weight = fh(0, 0.5)
  )

  # Convert ex1_delayed_effect to tte_data class
  ex1_delayed_effect_tte_data <- as.data.frame(ex1_delayed_effect)
  colnames(ex1_delayed_effect_tte_data) <- c("id", "tte", "event", "treatment")
  ex1_delayed_effect_tte_data$stratum <- "All"
  class(ex1_delayed_effect_tte_data) <- c("tte_data", "data.frame")

  expected <- wlr(ex1_delayed_effect_tte_data, weight = fh(0, 0.5))

  expect_equal(observed, expected)
})

test_that("wlr() accepts formula for stratified design", {
  data("ex1_delayed_effect", package = "simtrial", envir = environment())
  ex1_delayed_effect$trt <- ifelse(
    ex1_delayed_effect$trt == 1,
    "experimental",
    "control"
  )
  ex1_delayed_effect$strtm <- sample(
    x = c("biomarker positive", "biomarker negative"),
    size = nrow(ex1_delayed_effect),
    replace = TRUE,
    prob = c(0.6, 0.4)
  )

  # Use ex1_delayed_effect directly via the argument `formula`
  observed <- wlr(
    data = ex1_delayed_effect,
    formula = Surv(month, evntd) ~ trt + strata(strtm),
    weight = fh(0, 0.5)
  )

  # Convert ex1_delayed_effect to tte_data class
  ex1_delayed_effect_tte_data <- as.data.frame(ex1_delayed_effect)
  colnames(ex1_delayed_effect_tte_data) <- c("id", "tte", "event", "treatment", "stratum")
  class(ex1_delayed_effect_tte_data) <- c("tte_data", "data.frame")

  expected <- wlr(ex1_delayed_effect_tte_data, weight = fh(0, 0.5))

  expect_equal(observed, expected)
})

test_that("wlr() warns when formula argument is ignored", {
  x <- sim_pw_surv(n = 300) |> cut_data_by_event(100)
  expect_warning(
    wlr(x, weight = fh(0, 0.5), formula = Surv(tte, event) ~ treatment),
    "The formula argument was ignored"
  )

  y <- counting_process(x, arm = "experimental")
  expect_warning(
    wlr(y, weight = fh(0, 0.5), formula = Surv(tte, event) ~ treatment),
    "The formula argument was ignored"
  )
})

test_that("wlr.default() and wlr.tte_data() require arm='experimental'", {
  x <- sim_pw_surv(n = 300) |> cut_data_by_event(100)
  x$treatment <- ifelse(x$treatment == "experimental", "test", x$treatment)

  expect_error(
    wlr(x, weight = fh(0, 0.5)),
    "counting_process: arm is not a valid treatment group value."
  )
  expect_error(
    wlr(as.data.frame(x), weight = fh(0, 0.5), formula = Surv(tte, event) ~ treatment),
    "counting_process: arm is not a valid treatment group value."
  )

  # To use a different value for arm, have to manually run counting_process()
  y <- counting_process(x, arm = "test")
  expect_silent(wlr(y, weight = fh(0, 0.5)))
})
