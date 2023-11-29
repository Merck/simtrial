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
  two_na_var <- c(NA, NA)
  expect_error(
    input_check_scalar(two_na_var),
    "two_na_var must be a single non-negative number \\(or NA\\)"
  )
  two_num_var <- c(1, 2)
  expect_error(
    input_check_scalar(two_num_var),
    "two_num_var must be a single non-negative number \\(or NA\\)"
  )
  na_num_var <- c(NA, 2)
  expect_error(
    input_check_scalar(na_num_var),
    "na_num_var must be a single non-negative number \\(or NA\\)"
  )
})

test_that("input_check_scalar() fails non-numeric input", {
  logical_var <- TRUE
  expect_error(
    input_check_scalar(logical_var),
    "logical_var must be a single non-negative number \\(or NA\\)"
  )
  string_var <- "string"
  expect_error(
    input_check_scalar(string_var),
    "string_var must be a single non-negative number \\(or NA\\)"
  )
})

test_that("input_check_scalar() fails negative numbers", {
  neg_num_var <- -1
  expect_error(
    input_check_scalar(neg_num_var),
    "neg_num_var must be a single non-negative number \\(or NA\\)"
  )
  neg_num_vec_var <- c(-1, -2)
  expect_error(
    input_check_scalar(neg_num_vec_var),
    "neg_num_vec_var must be a single non-negative number \\(or NA\\)"
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
    "notwhole must be a single non-negative whole number \\(or NA\\)"
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
  with_zero_var <- 0:3
  expect_error(
    input_check_vector(with_zero_var),
    "with_zero_var must be a vector with only positive numbers and missing values"
  )
})

test_that("input_check_vector() fails with negative numbers", {
  negative_var <- -3:3
  expect_error(
    input_check_vector(negative_var),
    "negative_var must be a vector with only positive numbers and missing values"
  )
})

test_that("input_check_vector() fails with non-numeric vectors", {
  logical_var <- TRUE
  expect_error(
    input_check_vector(logical_var),
    "logical_var must be a vector with only positive numbers and missing values"
  )
  string_var <- "string"
  expect_error(
    input_check_vector(string_var),
    "string_var must be a vector with only positive numbers and missing values"
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
    "notwhole must be a vector with only positive whole numbers and missing values"
  )
})

# is_whole_number() ------------------------------------------------------------

test_that("is_whole_number() can distinguish between whole and decimal numbers", {
  x <- c(1.1, -1.1, 0, 2, NA)
  observed <- is_whole_number(x)
  expected <- c(FALSE, FALSE, TRUE, TRUE, NA)
  expect_equal(observed, expected)
})

test_that("is_whole_number() can tolerate minute computational inaccuracies from arithmetic operations", {
  expect_false(is_whole_number(3.0000001))
  expect_true(is_whole_number(3.00000001))

  expect_false(is_whole_number(3.000000000000001, tol = 1e-16))
  expect_true(is_whole_number(3.0000000000000001, tol = 1e-16))
})
