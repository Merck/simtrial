### bivariate with small var####
set.seed(1)
pbvnorm <- function(h1, h2, ro) {
  
  x = c(0.04691008, 0.23076534, 0.5, 0.76923466, 0.95308992)
  w = c(0.018854042, 0.038088059, 0.0452707394, 0.038088059, 0.018854042)
  h12 = (h1 * h1 + h2 * h2)/2
  bv = 0
  if (abs(ro) < 0.7) {
    h3 = h1 * h2
    ror = ro * x
    ror2 = 1 - ror * ror
    bv = ro * w * exp((ror * h3 - h12)/ror2)/sqrt(ror2)
    bvsum = sum(bv)
    bvfinal = bvsum + pnorm(h1, lower.tail = TRUE) * pnorm(h2, lower.tail = TRUE)
    
  } else {
    r2 = 1 - ro * ro
    r3 = sqrt(r2)
    if (ro <= 0) 
      h2 = -h2
    h3 = h1 * h2
    h7 = exp(-h3/2)
    if (r2 == 0) {
      if (ro > 0) 
        bvfinal = pnorm(min(h1, h2)) + bv * r3 * h7
      if (ro <= 0) 
        bvfinal = max(0, pnorm(h1) - pnorm(h2)) - bv * r3 * h7
    } else {
      
      h6 = abs(h1 - h2)
      h5 = h6 * h6/2
      h6 = h6/r3
      aa = 0.5 - h3/8
      ab = 3 - 2 * aa * h5
      bv = 0.13298076 * h6 * ab * pnorm(-h6) - exp(-h5/r2) * (ab + 
                                                                aa * r2) * 0.053051647
      r1 = r3 * x
      rr = r1 * r1
      r2 = sqrt(1 - rr)
      bv1 = w * exp(-h5/rr) * (exp(-h3/(1 + r2))/r2/h7 - 1 - aa * 
                                 rr)
      bvsum1 = sum(bv1)
      bv = bv - bvsum1
      if (ro > 0) 
        bvfinal = pnorm(min(h1, h2)) + bv * r3 * h7
      if (ro <= 0) 
        bvfinal = max(0, pnorm(h1) - pnorm(h2)) - bv * r3 * h7
    }
  }
}

testthat::test_that("bivariate with small var", {
  p1 = pbvnorm(1, 1, 0.5)
  rho = matrix(c(1, 0.5, 0.5, 1), 2, 2)
  p2 = mvtnorm::pmvnorm(-Inf, c(1, 1), c(0, 0), rho, algorithm = mvtnorm::GenzBretz(maxpts = 50000, 
                                                                           abseps = 1e-05))
  testthat::expect_equal(p1, p2[1], tolerance = 1e-05)
})

### bivariate with large var####
pbvnorm2 <- function(h1, h2, ro) {
  
  x = c(0.04691008, 0.23076534, 0.5, 0.76923466, 0.95308992)
  w = c(0.018854042, 0.038088059, 0.0452707394, 0.038088059, 0.018854042)
  h12 = (h1 * h1 + h2 * h2)/2
  bv = 0
  r2 = 1 - ro * ro
  r3 = sqrt(r2)
  h3 = h1 * h2
  h7 = exp(-h3/2)
  h6 = abs(h1 - h2)
  h5 = h6 * h6/2
  h6 = h6/r3
  aa = 0.5 - h3/8
  ab = 3 - 2 * aa * h5
  bv = 0.13298076 * h6 * ab * pnorm(h6) - exp(-h5/r2) * (ab + aa * r2) * 
    0.053051647
  r1 = r3 * x
  rr = r1 * r1
  r2 = sqrt(1 - rr)
  bv1 = w * exp(-h5/rr) * (exp(-h3/(1 + r2))/r2/h7 - 1 - aa * rr)
  bvsum1 = sum(bv1)
  bv = bv - bvsum1
  if (ro > 0) 
    bvfinal = pnorm(max(h1, h2)) + bv * r3 * h7
  if (ro <= 0) 
    bvfinal = max(0, pnorm(h1) - pnorm(h2)) - bv * r3 * h7
  bvfinal
}

testthat::test_that("bivariate with large var", {
  p1 = pbvnorm2(1, 1, 0.9)
  rho = matrix(c(1, 0.9, 0.9, 1), 2, 2)
  p2 = mvtnorm::pmvnorm(-Inf, c(1, 1), c(0, 0), rho, algorithm = mvtnorm::GenzBretz(maxpts = 50000, 
                                                                           abseps = 1e-05))
  testthat::expect_equal(p1, p2[1], tolerance = 1e-05)
})


### trivariate####
ptvnorm <- function(h, r, ro) {
  x = c(0.0491008, 0.23076534, 0.5, 0.76923466, 0.95308992)
  w = c(0.018854042, 0.038088059, 0.0452707394, 0.038088059, 0.018854042)
  
  h1 = h[1]
  h2 = h[2]
  h3 = h[3]
  r12 = r[1]
  r13 = r[2]
  r23 = r[3]
  
  if ((abs(r23) < abs(r12)) || (abs(r23) < abs(r13))) {
    if (abs(r12) >= abs(r13)) {
      hh = h1
      h1 = h3
      h3 = h2
      h2 = hh
      rr = r23
      r23 = r12
      r12 = r13
      r13 = rr
    } else {
      hh = h1
      h1 = h2
      h2 = hh
      rr = r23
      r23 = r13
      r13 = rr
    }
  }
  rh = c(h1, h2, h3, r12, r13, r23)
  h12 = h1 * h2
  h13 = h1 * h3
  h122 = (h1 * h1 + h2 * h2)/2
  h132 = (h1 * h1 + h3 * h3)/2
  del1 = 1 - r12 * r12 - r13 * r13 - r23 * r23 + 2 * r12 * r13 * r23
  
  rr12 = r12 * x
  rr13 = r13 * x
  del = 1 - rr12 * rr12 - rr13 * rr13 - r23 * r23 + 2 * rr12 * rr13 * 
    r23
  fac = sqrt(del)
  rr122 = 1 - rr12 * rr12
  rr133 = 1 - rr13 * rr13
  f1 = rr13 - r23 * rr12
  f2 = r23 - rr12 * rr13
  f3 = rr12 - r23 * rr13
  hp1 = (h3 * rr122 - h1 * f1 - h2 * f2)/fac/sqrt(rr122)
  hp2 = (h2 * rr133 - h1 * f3 - h3 * f2)/fac/sqrt(rr133)
  TV = w * exp((rr12 * h12 - h122)/rr122)/sqrt(rr122) * pnorm(hp1) * 
    r12 + w * exp((rr13 * h13 - h132)/rr133)/sqrt(rr133) * pnorm(hp2) * 
    r13
  TV = sum(TV)
  rho = matrix(c(1, r23, r23, 1), 2, 2)
  p2 = mvtnorm::pmvnorm(-Inf, c(h2, h3), c(0, 0), rho)
  
  TV = (TV + pnorm(h1) * p2)
  TV
}

testthat::test_that("trivariate", {
  corr1 = sqrt(344/550)
  corr2 = sqrt(190/550)
  corr3 = sqrt(139^2/(344 * 190))
  rho = matrix(c(1, corr1, corr2, corr1, 1, corr3, corr2, corr3, 1), 
               3, 3)
  p1 = ptvnorm(c(1, 1, 1), c(corr1, corr2, corr3))
  p2 = mvtnorm::pmvnorm(-Inf, c(1, 1, 1), c(0, 0, 0), rho, algorithm = mvtnorm::GenzBretz(maxpts = 50000, 
                                                                                 abseps = 1e-05))
  testthat::expect_equal(p1[1], p2[1], tolerance = 1e-05)
})

corr1 = sqrt(344/550)
corr2 = sqrt(190/550)
corr3 = sqrt(139^2/(344 * 190))
rho = matrix(c(1, corr1, corr2, corr1, 1, corr3, corr2, corr3, 1), 3, 3)
upper = c(-1, -1, -1)
ptvnorm(upper, c(corr1, corr2, corr3))[1]
mvtnorm::pmvnorm(-Inf, upper, c(0, 0, 0), rho, algorithm = Miwa())[1]
mvtnorm::pmvnorm(-Inf, upper, c(0, 0, 0), rho, algorithm = GenzBretz(maxpts = 50000, 
                                                                     abseps = 1e-07))[1]