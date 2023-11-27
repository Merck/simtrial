testthat::test_that("the z values match with the correspondings in fh_weight", {
  set.seed(1234)
  y <- sim_pw_surv(n = 300) |> cut_data_by_event(30)
  adjust.methods <- "asymp"
  wt <- list(a1 = c(0, 0), a2 = c(0, 1), a3 = c(1, 0), a4 = c(1, 1))
  ties.method <- "efron"
  one.sided <- TRUE
  HT.est <- FALSE
  max <- TRUE
  alpha <- 0.025
  data.anal <- data.frame(cbind(y$tte, y$event, y$treatment))
  fit <- survMisc::ten(survival::Surv(y$tte, y$event) ~ y$treatment, data = y)

  # Testing
  survMisc::comp(fit, p = sapply(wt, function(x) {
    x[1]
  }), q = sapply(wt, function(x) {
    x[2]
  }))
  tst.rslt <- attr(fit, "lrt")
  z1 <- tst.rslt$Z
  a2 <- y |> counting_process(arm = "experimental")
  aa <- fh_weight(a2, rho_gamma = data.frame(rho = c(0, 0, 1, 1), gamma = c(0, 1, 0, 1)))
  z2 <- aa$z
  expect_equal(c(z1[1], z1[7:9]), z2, tolerance = 0.00001)
})

testthat::test_that("fh_weight calculated correct correlation value when input a sequence of rho and gamma", {
  set.seed(123)
  y <- sim_pw_surv(n = 300) |> cut_data_by_event(30)
  adjust.methods <- "asymp"
  wt <- list(a1 = c(0, 0), a2 = c(0, 1), a3 = c(1, 0), a4 = c(1, 1))
  ties.method <- "efron"
  one.sided <- TRUE
  HT.est <- FALSE
  max <- TRUE
  alpha <- 0.025
  data.anal <- data.frame(cbind(y$tte, y$event, y$treatment))
  fit <- survMisc::ten(survival::Surv(y$tte, y$event) ~ y$treatment, data = y)

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
  fit2 <- survMisc::ten(survival::Surv(y$tte, y$event) ~ y$treatment, data = y)

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
  d1d2 <- dplyr::full_join(d1, d2, by = c("X1", "X2"))
  cov.tst.rslt1 <- d1d2$V
  cov.tst.1 <- matrix(NA, nrow = length(wt1), ncol = length(wt1))
  cov.tst.1[lower.tri(cov.tst.1, diag = FALSE)] <- cov.tst.rslt1
  cov.tst <- t(cov.tst.1)
  cov.tst[lower.tri(cov.tst, diag = FALSE)] <- cov.tst.rslt1
  diag(cov.tst) <- var.tst.rslt1
  cov.tst.1 <- matrix(Matrix::nearPD(cov.tst)$mat, length(Z.tst.rslt1), length(Z.tst.rslt1))
  z.max <- max(abs(tst.rslt1$Z))
  cor.tst <- cov2cor(cov.tst.1)
  pval2 <- 1 - max(min(mvtnorm::pmvnorm(
    lower = rep(-z.max, length(Z.tst.rslt1)),
    upper = rep(z.max, length(Z.tst.rslt1)),
    corr = cor.tst, algorithm = mvtnorm::Miwa()
  )[1], 0.9999), 0.0001)
  max.tst <- which(abs(Z.tst.rslt1) == max(abs(Z.tst.rslt1)), arr.ind = TRUE)
  if (Z.tst.rslt1[max.tst] >= 0) {
    pval <- 1 - pval2 / 2
  }
  if (Z.tst.rslt1[max.tst] < 0) {
    pval <- pval2 / 2
  }
  corr1 <- cor.tst[2:5, 2:5]
  a2 <- y |> counting_process(arm = "experimental")
  corr2 <- fh_weight(a2, rho_gamma = data.frame(rho = c(0, 0, 1, 1), gamma = c(0, 1, 0, 1)), return_corr = TRUE)
  corr2 <- rbind(corr2$v1, corr2$v2, corr2$v3, corr2$v4)
  expect_equal(corr1, corr2, tolerance = 0.00001)
})
