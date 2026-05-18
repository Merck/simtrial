testthat::test_that("tenFHcorr matches nphRCT Fleming-Harrington correlations",{
  testthat::skip_if_not_installed("nphRCT")

  set.seed(123)
  y <- sim_pw_surv(n = 300) %>% cut_data_by_event(30)
  rg <- tibble(rho = c(0, 0, 1, 1), gamma = c(0, 1, 0, 1))

  expected <- nphrct_tenFHcorr(y, rg)
  observed <- y %>%
    counting_process(arm = "Experimental") %>%
    tenFHcorr(rg = rg)

  expect_equal(observed$Z, expected$Z, tolerance = 0.00001)
  expect_equal(
    observed %>%
      dplyr::select(dplyr::starts_with("V")) %>%
      as.matrix(),
    expected %>%
      dplyr::select(dplyr::starts_with("V")) %>%
      as.matrix(),
    tolerance = 0.00001
  )
})
