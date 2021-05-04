
test_that("getCutDateForCount returns the correct cut date", {
  count = 100
  y = simPWSurv(n=200) 
  ycutdate = getCutDateForCount(y, count)
  
  x <- y %>% ungroup()%>% filter(fail==1) %>% arrange(cte) %>% slice(count)
  expect_equal(x$cte, ycutdate)
})