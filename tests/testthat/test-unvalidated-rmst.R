test_that("rmst() snapshot test", {
  data(ex1_delayed_effect)
  observed <- rmst(
    ex1_delayed_effect,
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
  data(ex1_delayed_effect)

  rmst_default <- rmst(
    ex1_delayed_effect,
    var_label_tte = "month",
    var_label_event = "evntd",
    var_label_group = "trt",
    tau = 10,
    reference = "0"
  )

  rmst_formula_1 <- rmst(
    month ~ evntd + trt,
    data = ex1_delayed_effect,
    tau = 10,
    reference = "0"
  )
  expect_equal(rmst_formula_1, rmst_default)

  rmst_formula_2 <- rmst(
    survival::Surv(month | evntd) ~ trt,
    data = ex1_delayed_effect,
    tau = 10,
    reference = "0"
  )
  expect_equal(rmst_formula_2, rmst_default)

  rmst_formula_3 <- rmst(
    ~ survival::Surv(month, evntd, trt),
    data = ex1_delayed_effect,
    tau = 10,
    reference = "0"
  )
  expect_equal(rmst_formula_3, rmst_default)
})
