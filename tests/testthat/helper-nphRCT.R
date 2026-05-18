nphrct_wlrt <- function(data, rho, gamma) {
  method <- if (rho == 0 && gamma == 0) {
    "lr"
  } else {
    "fh"
  }

  nphRCT::wlrt(
    formula = survival::Surv(tte, event) ~ Treatment,
    data = data,
    method = method,
    rho = if (method == "fh") rho else NULL,
    gamma = if (method == "fh") gamma else NULL,
    timefix = FALSE
  )
}

nphrct_wlr_z <- function(data, rg) {
  vapply(
    seq_len(nrow(rg)),
    function(i) nphrct_wlrt(data, rg$rho[i], rg$gamma[i])$z,
    numeric(1)
  )
}

nphrct_wlr_corr <- function(data, rg) {
  n_weight <- nrow(rg)
  cov_mat <- matrix(NA_real_, nrow = n_weight, ncol = n_weight)

  # FH covariance is the FH variance evaluated at averaged rho/gamma weights.
  for (i in seq_len(n_weight)) {
    for (j in seq_len(n_weight)) {
      cov_mat[i, j] <- nphrct_wlrt(
        data,
        rho = (rg$rho[i] + rg$rho[j]) / 2,
        gamma = (rg$gamma[i] + rg$gamma[j]) / 2
      )$v_u
    }
  }

  stats::cov2cor(cov_mat)
}

nphrct_tenFHcorr <- function(data, rg) {
  corr <- nphrct_wlr_corr(data, rg)
  colnames(corr) <- paste0("V", seq_len(ncol(corr)))
  tibble::as_tibble(cbind(rg, Z = nphrct_wlr_z(data, rg), corr))
}
