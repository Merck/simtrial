# 2024-02-23: Converted `example("maxcombo")` to tests from commit 914c604
# https://github.com/Merck/simtrial/pull/201/commits/914c6049cf9afb526f2286fdacfabf883050e175

test_that("maxcombo returns consistent results", {
  set.seed(1)
  observed <- sim_pw_surv(n = 200) |>
   cut_data_by_event(150) |>
   maxcombo(test1 = wlr(data, rho = 0, gamma = 0) |> quote(),
            test2 = wlr(data, rho = 0, gamma = 0.5) |> quote())
  expected <- data.frame(p_value = 1.5739680815363144e-06)
  expect_equal(observed, expected)
})
