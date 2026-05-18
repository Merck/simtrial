testthat::test_that("wlr Z values match nphRCT Fleming-Harrington tests",{
  testthat::skip_if_not_installed("nphRCT")

  set.seed(1234)
  y <- sim_pw_surv(n = 300) %>% cut_data_by_event(30)
  rg <- tibble(rho = c(0, 0, 1, 1), gamma = c(0, 1, 0, 1))

  expected_z <- nphrct_wlr_z(y, rg)
  observed_z <- y %>%
    counting_process(arm = "Experimental") %>%
    wlr(rg = rg) %>%
    dplyr::pull(Z)

  expect_equal(observed_z, expected_z, tolerance = 0.00001)
})
