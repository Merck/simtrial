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

test_that("input_check_scalar() can check for whole numbers", {
  whole <- 1
  missing <- NA
  notwhole <- 1.1

  expect_silent(input_check_scalar(whole))
  expect_silent(input_check_scalar(missing))
  expect_silent(input_check_scalar(notwhole))
  expect_silent(input_check_scalar(whole, require_whole_number = TRUE))
  expect_silent(input_check_scalar(missing, require_whole_number = TRUE))
  expect_error(
    input_check_scalar(notwhole, require_whole_number = TRUE),
    "x must be a single non-negative whole number \\(or NA\\)"
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

test_that("input_check_vector() can check for whole numbers", {
  whole <- c(1, 2, NA, 4, 5)
  notwhole <- c(1, 2.2, NA, 4, 5)

  expect_silent(input_check_vector(whole))
  expect_silent(input_check_vector(notwhole))
  expect_silent(input_check_vector(whole, require_whole_number = TRUE))
  expect_error(
    input_check_vector(notwhole, require_whole_number = TRUE),
    "x must be a vector with only positive whole numbers and missing values"
  )
})

# is_whole_number() ------------------------------------------------------------

test_that("is_whole_number() can distinguish between whole and decimal numbers", {
  x <- c(1.1, -1.1, 0, 2, NA)
  observed <- is_whole_number(x)
  expected <- c(FALSE, FALSE, TRUE, TRUE, NA)
  expect_equal(observed, expected)
})
