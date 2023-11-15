# input_check_scalar() ---------------------------------------------------------

test_that("input_check_scalar() passes a single non-negative number", {
  expect_silent(input_check_scalar(0))
  expect_silent(input_check_scalar(1))
  expect_silent(input_check_scalar(1L))
  expect_silent(input_check_scalar(100))
})

test_that("input_check_scalar() passes a single missing value", {
  expect_silent(input_check_scalar(NA))
})

test_that("input_check_scalar() fails multiple values", {
  expect_error(
    input_check_scalar(c(NA, NA)),
    "x must be a single non-negative number \\(or NA\\)"
  )
  expect_error(
    input_check_scalar(c(1, 2)),
    "x must be a single non-negative number \\(or NA\\)"
  )
  expect_error(
    input_check_scalar(c(NA, 2)),
    "x must be a single non-negative number \\(or NA\\)"
  )
})

test_that("input_check_scalar() fails non-numeric input", {
  expect_error(
    input_check_scalar(TRUE),
    "x must be a single non-negative number \\(or NA\\)"
  )
  expect_error(
    input_check_scalar("string"),
    "x must be a single non-negative number \\(or NA\\)"
  )
})

test_that("input_check_scalar() fails negative numbers", {
  expect_error(
    input_check_scalar(-1),
    "x must be a single non-negative number \\(or NA\\)"
  )
  expect_error(
    input_check_scalar(c(-1, -2)),
    "x must be a single non-negative number \\(or NA\\)"
  )
})

test_that("input_check_scalar() includes label in error message", {
  expect_error(
    input_check_scalar(-1, label = "custom"),
    "custom must be a single non-negative number \\(or NA\\)"
  )
})

# input_check_vector() ---------------------------------------------------------

test_that("input_check_vector() passes a vector of positive numbers", {
  expect_silent(input_check_vector(1:3))
  expect_silent(input_check_vector(c(1.1, 2.2, 3.3)))
})

test_that("input_check_vector() passes a vector of positive numbers and missing values", {
  expect_silent(input_check_vector(c(NA, 1)))
  expect_silent(input_check_vector(c(1, NA)))
})

test_that("input_check_vector() passes a single missing value", {
  expect_silent(input_check_vector(NA))
})

test_that("input_check_vector() passes a vector of missing values", {
  expect_silent(input_check_vector(c(NA, NA)))
})

test_that("input_check_vector() fails with zero", {
  expect_error(
    input_check_vector(0:3),
    "x must be a vector with only positive numbers and missing values"
  )
})

test_that("input_check_vector() fails with negative numbers", {
  expect_error(
    input_check_vector(-3:3),
    "x must be a vector with only positive numbers and missing values"
  )
})

test_that("input_check_vector() fails with non-numeric vectors", {
  expect_error(
    input_check_vector(TRUE),
    "x must be a vector with only positive numbers and missing values"
  )
  expect_error(
    input_check_vector("string"),
    "x must be a vector with only positive numbers and missing values"
  )
})

test_that("input_check_vector() includes label in error message", {
  expect_error(
    input_check_vector(-1, label = "custom"),
    "custom must be a vector with only positive numbers and missing values"
  )
})
