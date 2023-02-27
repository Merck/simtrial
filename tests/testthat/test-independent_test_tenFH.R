testthat::test_that("the Z values match with the correspondings in tenFH",{
  set.seed(1234)
  y=simPWSurv(n=300) %>% cutDataAtCount(30)
  adjust.methods="asymp"
  wt=list(a1=c(0,0),a2=c(0,1),a3=c(1,0),a4=c(1,1))
  ties.method = "efron"
  one.sided   = TRUE
  HT.est      = FALSE
  max         = TRUE
  alpha       = 0.025
  data.anal <- data.frame(cbind(y$tte,y$event,y$Treatment))
  fit<- survMisc::ten(Surv(y$tte, y$event) ~ y$Treatment, data = y)

  #Testing
  survMisc::comp(fit, p= sapply(wt, function(x){x[1]}), q= sapply(wt, function(x){x[2]}))
  tst.rslt <- attr(fit, 'lrt')
  z1=tst.rslt$Z
  a2 <- y %>% counting_process(arm="Experimental")
  aa=tenFH(a2,rg=tibble(rho=c(0,0,1,1),gamma=c(0,1,0,1)))
  z2=aa$Z
  expect_equal(c(z1[1],z1[7:9]),z2,tolerance = 0.00001)
})
