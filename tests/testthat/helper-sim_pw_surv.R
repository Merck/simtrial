# Helper functions used by test-double_programming_sim_pw_surv.R

test_sim_pw_surv <- function() {
  stratum <- data.frame(stratum = c("Low", "High"), p = c(.4, .6))

  block <- c(rep("control", 2), rep("experimental", 2))

  enroll_rate <- data.frame(duration = c(5, 195), rate = c(100, 3000))

  fail_rate <- dplyr::bind_rows(
    data.frame(stratum = "Low", period = 1, treatment = "control", duration = 3, rate = .03),
    data.frame(stratum = "Low", period = 2, treatment = "control", duration = 297, rate = .03),
    data.frame(stratum = "Low", period = 1, treatment = "experimental", duration = 3, rate = .03),
    data.frame(stratum = "Low", period = 2, treatment = "experimental", duration = 297, rate = .02),
    data.frame(stratum = "High", period = 1, treatment = "control", duration = 3, rate = .05),
    data.frame(stratum = "High", period = 2, treatment = "control", duration = 297, rate = .05),
    data.frame(stratum = "High", period = 1, treatment = "experimental", duration = 3, rate = .06),
    data.frame(stratum = "High", period = 2, treatment = "experimental", duration = 297, rate = .03)
  )
  dropout_rate <- dplyr::bind_rows(
    data.frame(stratum = "Low", period = 1, treatment = "control", duration = 300, rate = .001),
    data.frame(stratum = "Low", period = 1, treatment = "experimental", duration = 300, rate = .001),
    data.frame(stratum = "High", period = 1, treatment = "control", duration = 300, rate = .001),
    data.frame(stratum = "High", period = 1, treatment = "experimental", duration = 300, rate = .001)
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

  # Prepare to test block
  block1 <- x |> dplyr::filter(stratum == "Low")
  block2 <- x |> dplyr::filter(stratum == "High")
  bktest1 <- c()
  j <- 1
  for (i in seq(1, floor(nrow(block1) / 4))) {
    j <- 4 * i - 3
    bktest1[i] <- sum(str_count(block1$treatment[j:(j + 3)], "control"))
  }
  j <- 1
  bktest2 <- 0
  for (i in seq(1, floor(nrow(block2) / 4))) {
    j <- 4 * i - 3
    bktest2[i] <- sum(str_count(block2$treatment[j:(j + 3)], "control"))
  }

  # Prepare to test fail_rate
  y <- cut_data_by_date(x, cut_date = 300)

  intervals <- c(3)
  rate00 <- with(subset(y, treatment == "control" | stratum == "Low"), fit_pwexp(Surv(tte, event), intervals))
  rate01 <- with(subset(y, treatment == "control" | stratum == "High"), fit_pwexp(Surv(tte, event), intervals))
  rate10 <- with(subset(y, treatment == "experimental" | stratum == "Low"), fit_pwexp(Surv(tte, event), intervals))
  rate11 <- with(subset(y, treatment == "experimental" | stratum == "High"), fit_pwexp(Surv(tte, event), intervals))
  ratetest <- c(rate00$rate, rate10$rate, rate01$rate, rate11$rate)
  xevent <- dplyr::bind_rows(rate00, rate01, rate10, rate11)

  list(
    "x" = x,
    "bktest1" = bktest1,
    "bktest2" = bktest2,
    "ratetest" = ratetest,
    "stratum" = stratum,
    "block" = block,
    "enroll_rate" = enroll_rate,
    "fail_rate" = fail_rate,
    "dropout_rate" = dropout_rate,
    "xevent" = xevent
  )
}

# Check the arguments, by changing n, the actual number of events changes
test_sim_pw_surv_2 <- function() {
  res <- test_sim_pw_surv()

  set.seed(2468)
  z <- sim_pw_surv(
    n = 300000,
    stratum = res$stratum,
    block = res$block,
    enroll_rate = res$enroll_rate,
    fail_rate = res$fail_rate,
    dropout_rate = res$dropout_rate
  )

  y1 <- cut_data_by_date(z, cut_date = 300)

  intervals <- c(3)
  rate00 <- with(subset(y1, treatment == "control" | stratum == "Low"), fit_pwexp(Surv(tte, event), intervals))
  rate01 <- with(subset(y1, treatment == "control" | stratum == "High"), fit_pwexp(Surv(tte, event), intervals))
  rate10 <- with(subset(y1, treatment == "experimental" | stratum == "Low"), fit_pwexp(Surv(tte, event), intervals))
  rate11 <- with(subset(y1, treatment == "experimental" | stratum == "High"), fit_pwexp(Surv(tte, event), intervals))
  zevent <- dplyr::bind_rows(rate00, rate01, rate10, rate11)

  list("zevent" = zevent)
}
