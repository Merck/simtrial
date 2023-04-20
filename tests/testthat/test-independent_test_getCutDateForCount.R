
test_that("get_cut_date_by_event returns the correct cut date", {
  count = 100
  y = sim_pw_surv(n=200)
  ycutdate = get_cut_date_by_event(y, count)

  x <- y %>% ungroup()%>% filter(fail==1) %>% arrange(cte) %>% slice(count)
  expect_equal(x$cte, ycutdate)
})
