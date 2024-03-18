test_that("rmst() snapshot test", {
  data("ex1_delayed_effect")
  observed <- rmst(
    data = ex1_delayed_effect,
    var_label_tte = "month",
    var_label_event = "evntd",
    var_label_group = "trt",
    tau = 10,
    reference = "0"
  )
  expected <- data.frame(
    rmst_arm1 = 6.495175253205431,
    rmst_arm0 = 5.630125973237457,
    rmst_diff = 0.8650492799679741,
    z = 2.2178796367487963
  )
  expect_equal(observed, expected)
})

test_that("formula method matches default method", {
  data("ex1_delayed_effect")

  rmst_default <- rmst(
    data = ex1_delayed_effect,
    var_label_tte = "month",
    var_label_event = "evntd",
    var_label_group = "trt",
    tau = 10,
    reference = "0"
  )

  rmst_formula_1 <- rmst(
    data = ex1_delayed_effect,
    formula = month ~ evntd + trt,
    tau = 10,
    reference = "0"
  )
  expect_equal(rmst_formula_1, rmst_default)

  rmst_formula_2 <- rmst(
    data = ex1_delayed_effect,
    formula = survival::Surv(month | evntd) ~ trt,
    tau = 10,
    reference = "0"
  )
  expect_equal(rmst_formula_2, rmst_default)

  rmst_formula_3 <- rmst(
    data = ex1_delayed_effect,
    formula = ~ survival::Surv(month, evntd, trt),
    tau = 10,
    reference = "0"
  )
  expect_equal(rmst_formula_3, rmst_default)
})

test_that("default and formula methods of rmst are pipeable", {
  data("ex1_delayed_effect")

  rmst_default <- ex1_delayed_effect |>
    rmst(
      var_label_tte = "month",
      var_label_event = "evntd",
      var_label_group = "trt",
      tau = 10,
      reference = "0"
    )

  rmst_formula_1 <- ex1_delayed_effect |>
    rmst(
      formula = month ~ evntd + trt,
      tau = 10,
      reference = "0"
    )

  expect_equal(rmst_formula_1, rmst_default)
})

test_that("formula argument throws error for bad input data", {
  data("ex1_delayed_effect")

  # formula with 2 variables
  expect_error(
    rmst_formula_1 <- rmst(
      data = ex1_delayed_effect,
      formula = month ~ evntd,
      tau = 10,
      reference = "0"
    ),
    "The formula interface requires exactly 3 variables specified"
  )

  # formula with 4 variables
  expect_error(
    rmst_formula_1 <- rmst(
      data = ex1_delayed_effect,
      formula = month ~ evntd + trt + id,
      tau = 10,
      reference = "0"
    ),
    "The formula interface requires exactly 3 variables specified"
  )

  # non-formula
  expect_error(
    rmst_formula_1 <- rmst(
      data = ex1_delayed_effect,
      formula = "month ~ evntd + trt + id",
      tau = 10,
      reference = "0"
    ),
    'inherits(formula, "formula") is not TRUE',
    fixed = TRUE
  )
})
