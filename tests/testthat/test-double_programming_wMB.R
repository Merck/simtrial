test_mb_weight <- function(x, delay = 4) {
  out <- NULL
  for (i in unique(x$stratum)) {
    outi <- x[x$stratum == i, ]
    outi.sort <- outi[order(outi$tte), ]
    locmaxt <- length(outi.sort$tte[outi.sort$tte <= delay]) # location of the maximum timepoint (tte) that is less or equal to the input 'delay'
    outi$mb_weight <- NA
    outi$mb_weight[1:locmaxt] <- 1 / outi$s[1:locmaxt]
    outi$mb_weight[(locmaxt + 1):nrow(outi)] <- outi$mb_weight[locmaxt]
    out <- rbind(out, outi)
  }

  return(out)
}

# Test 1: for the situation of single stratum ####

test_that("Validation passed for the situation of single stratum", {
  x <- sim_pw_surv(n = 200) %>%
    cut_data_by_event(125) %>%
    counting_process(arm = "experimental")

  out1 <- test_mb_weight(x, delay = 3)
  out1 <- data.frame(out1[order(out1$stratum, out1$tte), ])
  out2 <- mb_weight(x, delay = 3)
  out2 <- data.frame(out2[order(out2$stratum, out2$tte), ])
  testthat::expect_equal(out1, out2)
})

# Test 2: for the situation of multiple stratum ####

test_that("Validation passed for the situation of multiple stratum", {
  set.seed(1)
  x <- sim_pw_surv(
    n = 200,
    # 2 stratum,30% and 70% prevalence
    stratum = tibble::tibble(stratum = c("Low", "High"), p = c(.3, .7)),
    fail_rate = tibble::tibble(
      stratum = c(rep("Low", 4), rep("High", 4)),
      period = rep(1:2, 4),
      treatment = rep(c(rep("control", 2), rep("experimental", 2)), 2),
      duration = rep(c(3, 1), 4),
      rate = c(.03, .05, .03, .03, .05, .08, .07, .04)
    ),
    dropout_rate = tibble::tibble(
      stratum = c(rep("Low", 2), rep("High", 2)),
      period = rep(1, 4),
      treatment = rep(c("control", "experimental"), 2),
      duration = rep(1, 4),
      rate = rep(.001, 4)
    )
  ) %>%
    cut_data_by_event(125) %>%
    counting_process(arm = "experimental")

  out1 <- test_mb_weight(x, delay = 3)
  out1 <- data.frame(out1[order(out1$stratum, out1$tte), ])
  out2 <- mb_weight(x, delay = 3)
  out2 <- data.frame(out2[order(out2$stratum, out2$tte), ])
  testthat::expect_equal(out1, out2)
})

# Test 3: for the situation where a stratum has no records before delay ends ####

test_that("Validation passed for the situation of a stratum with no records", {
  set.seed(1)
  x <- sim_pw_surv(
    n = 200,
    # 2 stratum,30% and 70% prevalence
    stratum = tibble::tibble(stratum = c("Low", "High"), p = c(.3, .7)),
    fail_rate = tibble::tibble(
      stratum = c(rep("Low", 4), rep("High", 4)),
      period = rep(1:2, 4),
      treatment = rep(c(rep("control", 2), rep("experimental", 2)), 2),
      duration = rep(c(3, 1), 4),
      rate = c(.03, .05, .03, .03, .05, .08, .07, .04)
    ),
    dropout_rate = tibble::tibble(
      stratum = c(rep("Low", 2), rep("High", 2)),
      period = rep(1, 4),
      treatment = rep(c("control", "experimental"), 2),
      duration = rep(1, 4),
      rate = rep(.001, 4)
    )
  ) %>%
    cut_data_by_event(125) %>%
    counting_process(arm = "experimental")

  out <- mb_weight(x, delay = 1)
  testthat::expect_false(anyNA(out$mb_weight))
})
