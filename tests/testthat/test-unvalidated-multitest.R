test_that("multitest() is equivalent to running tests individually", {
  trial_data <- sim_pw_surv(n = 200)
  trial_data_cut <- cut_data_by_event(trial_data, 150)

  # create cutting test functions
  wlr_partial <- create_test(wlr, weight = fh(rho = 0, gamma = 0))
  rmst_partial <- create_test(rmst, tau = 20)
  maxcombo_partial <- create_test(maxcombo, rho = c(0, 0), gamma = c(0, 0.5))

  skip('skip')
  observed <- multitest(
    data = trial_data_cut,
    wlr = wlr_partial,
    rmst = rmst_partial,
    maxcombo = maxcombo_partial
  )
  expected <- list(
    wlr = wlr(trial_data_cut, weight = fh(rho = 0, gamma = 0)),
    rmst = rmst(trial_data_cut, tau = 20),
    maxcombo = maxcombo(trial_data_cut, rho = c(0, 0), gamma = c(0, 0.5))
  )
  expect_equal(observed, expected)
})
