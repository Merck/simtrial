testthat::test_that("the Z values match with the correspondings in wlr",{
  set.seed(1234)
  y=sim_pw_surv(n=300) %>% cut_data_by_event(30)
  adjust.methods="asymp"
  wt=list(a1=c(0,0),a2=c(0,1),a3=c(1,0),a4=c(1,1))
  ties.method = "efron"
  one.sided   = TRUE
  HT.est      = FALSE
  max         = TRUE
  alpha       = 0.025
  data.anal <- data.frame(cbind(y$tte,y$event,y$treatment))
  fit<- survMisc::ten(Surv(y$tte, y$event) ~ y$treatment, data = y)

  #Testing
  survMisc::comp(fit, p= sapply(wt, function(x){x[1]}), q= sapply(wt, function(x){x[2]}))
  tst.rslt <- attr(fit, 'lrt')
  z1=tst.rslt$Z
  a2 <- y %>% counting_process(arm="experimental")
  aa=wlr(a2,rho_gamma=tibble(rho=c(0,0,1,1),gamma=c(0,1,0,1)))
  z2=aa$Z
  expect_equal(c(z1[1],z1[7:9]),z2,tolerance = 0.00001)
})
