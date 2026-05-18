testthat::test_that("pvalue_maxcombo matches nphRCT-based FH correlations",{
  testthat::skip_if_not_installed("nphRCT")

  set.seed(2022)
  y <- sim_pw_surv(n = 300) %>% cut_data_by_event(30)
  rg <- tibble(rho = c(0, 0, 1, 1), gamma = c(0, 1, 0, 1))

  z_expected <- nphrct_tenFHcorr(y, rg)
  z_observed <- y %>%
    counting_process(arm = "Experimental") %>%
    tenFHcorr(rg = rg)

  set.seed(1)
  expected_p <- 1 - mvtnorm::pmvnorm(
    lower = rep(min(z_expected$Z), nrow(z_expected)),
    corr = z_expected %>%
      dplyr::select(dplyr::starts_with("V")) %>%
      data.matrix(),
    algorithm = mvtnorm::GenzBretz(maxpts = 50000, abseps = 0.00001)
  )[1]

  set.seed(1)
  observed_p <- pvalue_maxcombo(
    Z = z_observed,
    algorithm = mvtnorm::GenzBretz(maxpts = 50000, abseps = 0.00001)
  )

  expect_equal(observed_p, as.numeric(expected_p), tolerance = 0.00001)
})
