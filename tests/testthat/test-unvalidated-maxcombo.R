# 2024-02-23: Converted `example("maxcombo")` to tests from commit 914c604
# https://github.com/Merck/simtrial/pull/201/commits/914c6049cf9afb526f2286fdacfabf883050e175

test_that("maxcombo returns consistent results", {
  set.seed(1)
  observed <- sim_pw_surv(n = 200) |>
    cut_data_by_event(150) |>
    maxcombo(rho = c(0, 0), gamma = c(0, 0.5))

  expect_equal(observed$p_value, 1.5739680815363144e-06)
})

test_that("maxcombo fails early with bad input", {
  input <- observed <- sim_pw_surv(n = 200) |> cut_data_by_event(150)

  expect_error(maxcombo(input, rho = c(-1, 0), gamma = c(0, 0.5)))
  expect_error(maxcombo(input, rho = c(0, 0), gamma = c(-1, 0.5)))
  expect_error(maxcombo(input, rho = letters[1:2], gamma = c(0, 0.5)))
  expect_error(maxcombo(input, rho = c(0, 0), gamma = letters[1:2]))
  expect_error(maxcombo(input, rho = c(0), gamma = c(0, 0.5)))
})
