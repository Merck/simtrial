stratum <- tibble::tibble(stratum = c("Low", "High"), p = c(.4, .6))

block <- c(rep("control", 2), rep("experimental", 2))

enroll_rate <- tibble::tibble(duration = c(5, 195), rate = c(100, 3000))

fail_rate <- dplyr::bind_rows(
  tibble::tibble(stratum = "Low", period = 1, treatment = "control", duration = 3, rate = .03),
  tibble::tibble(stratum = "Low", period = 2, treatment = "control", duration = 297, rate = .03),
  tibble::tibble(stratum = "Low", period = 1, treatment = "experimental", duration = 3, rate = .03),
  tibble::tibble(stratum = "Low", period = 2, treatment = "experimental", duration = 297, rate = .02),
  tibble::tibble(stratum = "High", period = 1, treatment = "control", duration = 3, rate = .05),
  tibble::tibble(stratum = "High", period = 2, treatment = "control", duration = 297, rate = .05),
  tibble::tibble(stratum = "High", period = 1, treatment = "experimental", duration = 3, rate = .06),
  tibble::tibble(stratum = "High", period = 2, treatment = "experimental", duration = 297, rate = .03)
)
dropout_rate <- dplyr::bind_rows(
  tibble::tibble(stratum = "Low", period = 1, treatment = "control", duration = 300, rate = .001),
  tibble::tibble(stratum = "Low", period = 1, treatment = "experimental", duration = 300, rate = .001),
  tibble::tibble(stratum = "High", period = 1, treatment = "control", duration = 300, rate = .001),
  tibble::tibble(stratum = "High", period = 1, treatment = "experimental", duration = 300, rate = .001)
)
set.seed(1)
x <- sim_pw_surv(
  n = 400000,
  stratum = stratum,
  block = block,
  enroll_rate = enroll_rate,
  fail_rate = fail_rate,
  dropout_rate = dropout_rate
)

# prepare to test block
block1 <- x %>% dplyr::filter(stratum == "Low")
block2 <- x %>% dplyr::filter(stratum == "High")
bktest1 <- c()
j <- 1
for (i in seq(1, floor(nrow(block1) / 4))) {
  j <- 4 * i - 3
  bktest1[i] <- sum(stringr::str_count(block1$treatment[j:(j + 3)], "control"))
}
j <- 1
bktest2 <- 0
for (i in seq(1, floor(nrow(block2) / 4))) {
  j <- 4 * i - 3
  bktest2[i] <- sum(stringr::str_count(block2$treatment[j:(j + 3)], "control"))
}

# prepare to test fail_rate

y <- cut_data_by_date(x, cut_date = 300)

intervals <- c(3)
rate00 <- with(subset(y, treatment == "control" | stratum == "Low"), fit_pwexp(Surv(tte, event), intervals))
rate01 <- with(subset(y, treatment == "control" | stratum == "High"), fit_pwexp(Surv(tte, event), intervals))
rate10 <- with(subset(y, treatment == "experimental" | stratum == "Low"), fit_pwexp(Surv(tte, event), intervals))
rate11 <- with(subset(y, treatment == "experimental" | stratum == "High"), fit_pwexp(Surv(tte, event), intervals))
ratetest <- c(rate00$rate, rate10$rate, rate01$rate, rate11$rate)
xevent <- dplyr::bind_rows(rate00, rate01, rate10, rate11)

testthat::test_that("stratum percentage calculated from simulated dataset must be within
                    the tolerance=0.002 of stratum in setup (0.4,0.6)", {
  expect_equal(
    object = c(
      sum(stringr::str_count(x$stratum, "Low")) / 400000,
      sum(stringr::str_count(x$stratum, "High")) / 400000
    ),
    expected = c(0.4, 0.6), tolerance = 0.002
  )
})

testthat::test_that("block calculated from simulated dataset equals size of 4 with 1:1
                    randomization, which is 2 for each arm", {
  expect_equal(object = bktest1, expected = rep(2, length(bktest1)))
  expect_equal(object = bktest2, expected = rep(2, length(bktest2)))
})

testthat::test_that("fail_rate calculated from simulated dataset must be within the
                    tolerance=0.1 of fail_rate in setting", {
  expect_equal(object = ratetest, expected = fail_rate$rate, tolerance = 0.1)
})

testthat::test_that("dropout_rate calculated from simulated dataset must be within
                    the tolerance=0.0005 of dropout_rate=0.001 in setup", {
  duration <- 300
  drtest <- 0
  for (i in 1:duration) {
    drtest[i] <- sum(x$dropout_time <= i & x$dropout_time > (i - 1)) / 400000
  }
  expect_equal(object = drtest, expected = rep(0.001, 300), tolerance = 0.001)
})

testthat::test_that("enroll_rate calculated from simulated dataset must be within
                    the relative tolerance=0.05 of enroll_rate in setup", {
  duration <- 300
  entest <- 0
  for (i in 1:duration) {
    entest[i] <- sum(x$enroll_time <= i & x$enroll_time > (i - 1))
  }
  entest1 <- entest[entest != 0]
  entestexp <- c(rep(100, 5), rep(3000, length(entest1) - 5))
  entest2 <- (entest1 - entestexp) / entestexp
  expect_equal(object = entest2, expected = rep(0, length(entest1)), tolerance = 0.05)
})

# check the arguments, by changing n, the actual number of events changes
set.seed(2468)
z <- sim_pw_surv(
  n = 300000,
  stratum = stratum,
  block = block,
  enroll_rate = enroll_rate,
  fail_rate = fail_rate,
  dropout_rate = dropout_rate
)


y1 <- cut_data_by_date(z, cut_date = 300)

intervals <- c(3)
rate00 <- with(subset(y1, treatment == "control" | stratum == "Low"), fit_pwexp(Surv(tte, event), intervals))
rate01 <- with(subset(y1, treatment == "control" | stratum == "High"), fit_pwexp(Surv(tte, event), intervals))
rate10 <- with(subset(y1, treatment == "experimental" | stratum == "Low"), fit_pwexp(Surv(tte, event), intervals))
rate11 <- with(subset(y1, treatment == "experimental" | stratum == "High"), fit_pwexp(Surv(tte, event), intervals))
zevent <- dplyr::bind_rows(rate00, rate01, rate10, rate11)

testthat::test_that("The actual number of events changes by changing total sample size", {
  expect_false(unique(xevent$event == zevent$event))
})
