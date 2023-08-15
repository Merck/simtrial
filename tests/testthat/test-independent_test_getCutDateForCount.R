test_that("get_cut_date_by_event returns the correct cut date", {
  event <- 100
  y <- sim_pw_surv(n = 200)
  ycutdate <- get_cut_date_by_event(y, event)

  x <- y %>%
    dplyr::ungroup() %>%
    filter(fail == 1) %>%
    dplyr::arrange(cte) %>%
    dplyr::slice(event)
  expect_equal(x$cte, ycutdate)
})
