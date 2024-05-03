test_that("test for sim_fixed_n power comparing to gsDesign results with fixed duration in timing_type=1", {
  skip_if_not_installed("gsDesign")

  test2 <- test_simfix()$test2
  tt1test <- subset(test2, test2$cut == "Planned duration", select = c(event, ln_hr, z, duration, sim))
  expect_equal(object = sum(as.integer(tt1test$z < (-1.96))) / 100, expected = 0.94, tolerance = 0.02)
})

test_that("test for sim_fixed_n power comparing to gsDesign results with target events in timing_type=2", {
  skip_if_not_installed("gsDesign")

  test2 <- test_simfix()$test2
  tt2test <- subset(test2, test2$cut == "Targeted events", select = c(event, ln_hr, z, duration, sim))
  expect_equal(object = sum(as.integer(tt2test$z < (-1.96))) / 100, expected = 0.93, tolerance = 0.02)
})

test_that("test for events in the correct directions in timing_type=3 comparing to timing_type=2", {
  skip_if_not_installed("gsDesign")

  test2 <- test_simfix()$test2

  tt2test <- subset(test2, test2$cut == "Targeted events", select = c(event, ln_hr, z, duration, sim))
  tt3test <- subset(test2, test2$cut == "Minimum follow-up", select = c(event, ln_hr, z, duration, sim))
  ttvalue <- 0
  for (i in seq_len(nrow(tt3test))) {
    if ((tt3test$duration[i] > tt2test$duration[i]) & (tt3test$event[i] >= tt2test$event[i])) {
      ttvalue[i] <- 1
    } else if ((tt3test$duration[i] <= tt2test$duration[i]) & (tt3test$event[i] <= tt2test$event[i])) {
      ttvalue[i] <- 1
    } else {
      ttvalue[i] <- 0
    }
  }

  expect_equal(object = unique(ttvalue), expected = 1)
})

test_that("test for timing_type=4 outputs using timing_type 1 and 2 output", {
  skip_if_not_installed("gsDesign")

  test2 <- test_simfix()$test2

  tt1test <- subset(test2, test2$cut == "Planned duration", select = c(event, ln_hr, z, duration, sim))
  tt2test <- subset(test2, test2$cut == "Targeted events", select = c(event, ln_hr, z, duration, sim))
  tt4test <- subset(test2, test2$cut == "Max(planned duration, event cut)", select = c(event, ln_hr, z, duration, sim))
  tt4event <- 0
  for (i in seq_len(nrow(tt4test))) {
    if (tt1test$duration[i] < tt2test$duration[i]) {
      tt4event[i] <- tt2test$event[i]
    } else {
      tt4event[i] <- tt1test$event[i]
    }
  }

  expect_equal(object = tt4event, expected = tt4test$event)
})

test_that("test for timing_type=5 outputs using timing_type 2 and 3 output", {
  skip_if_not_installed("gsDesign")

  test2 <- test_simfix()$test2

  tt2test <- subset(test2, test2$cut == "Targeted events", select = c(event, ln_hr, z, duration, sim))
  tt3test <- subset(test2, test2$cut == "Minimum follow-up", select = c(event, ln_hr, z, duration, sim))
  tt5test <- subset(test2, test2$cut == "Max(min follow-up, event cut)", select = c(event, ln_hr, z, duration, sim))
  tt5event <- 0
  for (i in seq_len(nrow(tt5test))) {
    if (tt2test$duration[i] < tt3test$duration[i]) {
      tt5event[i] <- tt3test$event[i]
    } else {
      tt5event[i] <- tt2test$event[i]
    }
  }

  expect_equal(object = tt5event, expected = tt5test$event)
})
