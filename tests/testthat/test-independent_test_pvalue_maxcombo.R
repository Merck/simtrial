testthat::test_that("the p-values correspond to pvalue_maxcombo", {
  set.seed(2022)
  # this part is a double programming
  y <- sim_pw_surv(n = 300) %>% cut_data_by_event(30)
  adjust.methods <- "asymp"
  wt <- list(a1 = c(0, 0), a2 = c(0, 1), a3 = c(1, 0), a4 = c(1, 1))
  ties.method <- "efron"
  one.sided <- TRUE
  HT.est <- FALSE
  max <- TRUE
  alpha <- 0.025
  fit <- survMisc::ten(Surv(y$tte, y$event) ~ y$treatment, data = y)

  # Testing
  survMisc::comp(fit, p = sapply(wt, function(x) {
    x[1]
  }), q = sapply(wt, function(x) {
    x[2]
  }))
  tst.rslt <- attr(fit, "lrt")

  # Combination test ("asymp")
  # Calculating the covariace matrix
  tst.rslt1 <- rbind(tst.rslt[1, ], subset(tst.rslt, grepl("FH", tst.rslt$W)))
  Z.tst.rslt1 <- tst.rslt1$Z
  q.tst.rslt1 <- tst.rslt1$Q
  var.tst.rslt1 <- tst.rslt1$Var
  wt1 <- c(list(a0 = c(0, 0)), wt)
  combo.wt <- combn(wt1, 2)
  combo.wt.list <- list()
  for (i in 1:ncol(combo.wt)) {
    combo.wt.list[[i]] <- combo.wt[, i]
  }
  combo.wt.list.up <- lapply(combo.wt.list, function(a) {
    mapply("+", a)
  })
  wt2 <- lapply(combo.wt.list.up, function(a) {
    apply(a, 1, "sum") / 2
  })
  d1 <- data.frame(do.call(rbind, wt2))
  wt3 <- unique(wt2)
  d2 <- data.frame(do.call(rbind, wt3))
  fit2 <- survMisc::ten(Surv(y$tte, y$event) ~ y$treatment, data = y)

  # Testing (for calculating the covariances)
  survMisc::comp(fit2, p = sapply(wt3, function(x) {
    x[1]
  }), q = sapply(wt3, function(x) {
    x[2]
  }))
  tst.rsltt <- attr(fit2, "lrt")
  tst.rslt2 <- subset(tst.rsltt, grepl("FH", tst.rsltt$W))
  cov.tst.rslt11 <- tst.rslt2$Var
  d2$V <- cov.tst.rslt11
  d1d2 <- full_join(d1, d2, by = c("X1", "X2"))
  cov.tst.rslt1 <- d1d2$V
  cov.tst.1 <- matrix(NA, nrow = length(wt1), ncol = length(wt1))
  cov.tst.1[lower.tri(cov.tst.1, diag = FALSE)] <- cov.tst.rslt1
  cov.tst <- t(cov.tst.1)
  cov.tst[lower.tri(cov.tst, diag = FALSE)] <- cov.tst.rslt1
  diag(cov.tst) <- var.tst.rslt1
  cov.tst.1 <- matrix(Matrix::nearPD(cov.tst)$mat, length(Z.tst.rslt1), length(Z.tst.rslt1))
  z.max <- max(abs(tst.rslt1$Z))
  cor.tst <- cov2cor(cov.tst.1)
  pval2 <- 1 - max(min(pmvnorm(
    lower = rep(-z.max, length(Z.tst.rslt1)),
    upper = rep(z.max, length(Z.tst.rslt1)),
    corr = cor.tst, algorithm = GenzBretz(maxpts = 50000, abseps = 0.00001)
  )[1], 0.9999), 0.0001)
  max.tst <- which(abs(Z.tst.rslt1) == max(abs(Z.tst.rslt1)), arr.ind = TRUE)
  if (Z.tst.rslt1[max.tst] >= 0) {
    pval <- 1 - pval2 / 2
  }
  if (Z.tst.rslt1[max.tst] < 0) {
    pval <- pval2 / 2
  }
  p1 <- pval

  a2 <- y %>% counting_process(arm = "experimental")
  aa <- tenFHcorr(a2, rho_gamma = tibble(rho = c(0, 0, 1, 1), gamma = c(0, 1, 0, 1)))
  p2 <- simtrial::pvalue_maxcombo(z = aa)

  expect_equal(p1, p2, tolerance = 0.005)
})
