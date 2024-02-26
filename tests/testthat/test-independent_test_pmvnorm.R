test_that("bivariate with small variance", {
  set.seed(1)

  p1 <- pbvnorm(1, 1, 0.5)
  rho <- matrix(c(1, 0.5, 0.5, 1), 2, 2)
  p2 <- mvtnorm::pmvnorm(
    lower = -Inf,
    upper = c(1, 1),
    mean = c(0, 0),
    corr = rho,
    algorithm = mvtnorm::GenzBretz(
      maxpts = 50000,
      abseps = 1e-05
    )
  )
  expect_equal(p1, p2[1], tolerance = 1e-05)
})

test_that("bivariate with large variance", {
  set.seed(1)

  p1 <- pbvnorm2(1, 1, 0.9)
  rho <- matrix(c(1, 0.9, 0.9, 1), 2, 2)
  p2 <- mvtnorm::pmvnorm(
    lower = -Inf,
    upper = c(1, 1),
    mean = c(0, 0),
    corr = rho,
    algorithm = mvtnorm::GenzBretz(
      maxpts = 50000,
      abseps = 1e-05
    )
  )
  expect_equal(p1, p2[1], tolerance = 1e-05)
})

test_that("trivariate", {
  set.seed(1)

  corr1 <- sqrt(344 / 550)
  corr2 <- sqrt(190 / 550)
  corr3 <- sqrt(139^2 / (344 * 190))
  rho <- matrix(
    c(1, corr1, corr2, corr1, 1, corr3, corr2, corr3, 1),
    3, 3
  )
  p1 <- ptvnorm(c(1, 1, 1), c(corr1, corr2, corr3))
  p2 <- mvtnorm::pmvnorm(
    lower = -Inf,
    upper = c(1, 1, 1),
    mean = c(0, 0, 0),
    corr = rho,
    algorithm = mvtnorm::GenzBretz(
      maxpts = 50000,
      abseps = 1e-05
    )
  )
  expect_equal(p1[1], p2[1], tolerance = 1e-05)
})

# corr1 <- sqrt(344 / 550)
# corr2 <- sqrt(190 / 550)
# corr3 <- sqrt(139^2 / (344 * 190))
# rho <- matrix(c(1, corr1, corr2, corr1, 1, corr3, corr2, corr3, 1), 3, 3)
# upper <- c(-1, -1, -1)
# ptvnorm(upper, c(corr1, corr2, corr3))[1]
# mvtnorm::pmvnorm(-Inf, upper, c(0, 0, 0), rho, algorithm = mvtnorm::Miwa())[1]
# mvtnorm::pmvnorm(-Inf, upper, c(0, 0, 0), rho, algorithm = mvtnorm::GenzBretz(
#   maxpts = 50000,
#   abseps = 1e-07
# ))[1]
