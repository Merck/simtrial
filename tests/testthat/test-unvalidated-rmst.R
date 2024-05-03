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
  expected <- list(
    rmst_diff = 0.8650492799679741,
    z = 2.2178796367487963
  )
  expect_equal(observed$estimation, expected$rmst_diff)
  expect_equal(observed$z, expected$z)
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

  rmst_formula <- rmst(
    data = ex1_delayed_effect,
    formula = Surv(month, evntd) ~ trt,
    tau = 10,
    reference = "0"
  )
  expect_equal(rmst_formula, rmst_default)
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
      formula = Surv(month, evntd) ~ trt,
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
      formula = Surv(month) ~ trt,
      tau = 10,
      reference = "0"
    ),
    "The formula interface requires exactly 3 variables specified"
  )
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
      formula = Surv(month, evntd) ~ trt + id,
      tau = 10,
      reference = "0"
    ),
    "The formula interface requires exactly 3 variables specified"
  )
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

  # non-canonical formula
  expect_error(
    rmst_formula_1 <- rmst(
      data = ex1_delayed_effect,
      formula = Surv(month | evntd) ~ trt,
      tau = 10,
      reference = "0"
    ),
    "Unable to identify a single tte variable"
  )
  expect_error(
    rmst_formula_1 <- rmst(
      data = ex1_delayed_effect,
      formula = ~ Surv(month, evntd, trt),
      tau = 10,
      reference = "0"
    ),
    "unused argument (trt)",
    fixed = TRUE
  )
  expect_error(
    rmst_formula_1 <- rmst(
      data = ex1_delayed_effect,
      formula = month ~ evntd + trt,
      tau = 10,
      reference = "0"
    ),
    "Must use canonical formula syntax with Surv()"
  )
  expect_error(
    rmst_formula_1 <- rmst(
      data = ex1_delayed_effect,
      formula = ~ month + evntd + trt,
      tau = 10,
      reference = "0"
    ),
    "Must use canonical formula syntax with Surv()"
  )
})

test_that("parse_formula_rmst() properly parses the formula argument", {
  expected <- c(
    "var_label_tte" = "tte",
    "var_label_event" = "event",
    "var_label_group" = "group"
  )

  expect_identical(
    parse_formula_rmst(formula = Surv(tte, event) ~ group),
    expected
  )

  expect_identical(
    parse_formula_rmst(formula = Surv(event = event, time = tte) ~ group),
    expected
  )


  expect_identical(
    parse_formula_rmst(formula = Surv(tte, event = event) ~ group),
    expected
  )

  expect_identical(
    parse_formula_rmst(formula = Surv(event = event, tte) ~ group),
    expected
  )

  # Note: 4 variables is not currently allowed. This invalid formula would be
  # caught upstream in rmst(). This test is just to show that
  # parse_formula_rmst() can still parse it correctly regardless.
  expect_identical(
    parse_formula_rmst(formula = Surv(tte, event) ~ group + group2),
    expected
  )
})
