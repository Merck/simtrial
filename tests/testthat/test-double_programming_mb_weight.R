double_progamming_mb_weight <- function(x, delay = 4) {
  out <- NULL
  for (i in unique(x$stratum)) {
    outi <- x[x$stratum == i, ]
    outi.sort <- outi[order(outi$tte), ]
    # location of the maximum timepoint (tte) that is less or equal to the input 'delay'
    locmaxt <- length(outi.sort$tte[outi.sort$tte <= delay])
    outi$weight <- NA
    outi$weight[1:locmaxt] <- 1 / outi$s[1:locmaxt]
    outi$weight[(locmaxt + 1):nrow(outi)] <- outi$weight[locmaxt]
    out <- rbind(out, outi)
  }

  return(out)
}


test_that("mb_weight works for single stratum", {
  x <- sim_pw_surv(n = 200) |>
    cut_data_by_event(125) |>
    counting_process(arm = "experimental")

  out1 <- double_progamming_mb_weight(x, delay = 3)
  out1 <- data.frame(out1[order(out1$stratum, out1$tte), ])
  out2 <- simtrial:::mb_weight(x, delay = 3)
  out2 <- data.frame(out2[order(out2$stratum, out2$tte), ])
  expect_equal(out1, out2)
})

test_that("mb_weight works for multiple strata", {
  set.seed(1)
  x <- sim_pw_surv(
    n = 200,
    # 2 stratum,30% and 70% prevalence
    stratum = data.frame(stratum = c("Low", "High"), p = c(.3, .7)),
    fail_rate = data.frame(
      stratum = c(rep("Low", 4), rep("High", 4)),
      period = rep(1:2, 4),
      treatment = rep(c(rep("control", 2), rep("experimental", 2)), 2),
      duration = rep(c(3, 1), 4),
      rate = c(.03, .05, .03, .03, .05, .08, .07, .04)
    ),
    dropout_rate = data.frame(
      stratum = c(rep("Low", 2), rep("High", 2)),
      period = rep(1, 4),
      treatment = rep(c("control", "experimental"), 2),
      duration = rep(1, 4),
      rate = rep(.001, 4)
    )
  ) |>
    cut_data_by_event(125) |>
    counting_process(arm = "experimental")

  out1 <- double_progamming_mb_weight(x, delay = 3)
  out1 <- data.frame(out1[order(out1$stratum, out1$tte), ])
  out2 <- simtrial:::mb_weight(x, delay = 3)
  out2 <- data.frame(out2[order(out2$stratum, out2$tte), ])
  expect_equal(out1, out2)
})

test_that("mb_weight works for a stratum with no records before delay ends", {
  set.seed(1)
  x <- sim_pw_surv(
    n = 200,
    # 2 stratum,30% and 70% prevalence
    stratum = data.frame(stratum = c("Low", "High"), p = c(.3, .7)),
    fail_rate = data.frame(
      stratum = c(rep("Low", 4), rep("High", 4)),
      period = rep(1:2, 4),
      treatment = rep(c(rep("control", 2), rep("experimental", 2)), 2),
      duration = rep(c(3, 1), 4),
      rate = c(.03, .05, .03, .03, .05, .08, .07, .04)
    ),
    dropout_rate = data.frame(
      stratum = c(rep("Low", 2), rep("High", 2)),
      period = rep(1, 4),
      treatment = rep(c("control", "experimental"), 2),
      duration = rep(1, 4),
      rate = rep(.001, 4)
    )
  ) |>
    cut_data_by_event(125) |>
    counting_process(arm = "experimental")

  out <- simtrial:::mb_weight(x, delay = 1)
  expect_false(anyNA(out$mb_weight))
})
