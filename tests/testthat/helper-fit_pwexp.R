# Helper functions used by test-double_programming_fit_pwexp.R

test_fit_pwexp <- function(Srv, intervals) {
  time <- Surv[, "time"]
  status <- Surv[, "status"]
  if (tail(time, 1) > sum(intervals) && tail(status, 1) == 1) intervals <- c(intervals, Inf)
  out <- NULL
  interval.start <- 0
  interval.end <- 0
  for (i in seq_along(intervals)) {
    if (i == length(intervals)) {
      interval.start <- 0
    } else {
      interval.start <- interval.start + intervals[i - 1]
    }
    interval.end <- interval.end + intervals[i]
    datai <- Srv[Srv[, "time"] > interval.start]
    if (nrow(datai) == 0) next
    datai[datai[, "time"] > interval.end][, "time"] <- interval.end
    datai[, "time"] <- datai[, "time"] - interval.start
    events <- sum(datai[, "status"])
    TTOT <- sum(datai[, "times"])
    out <- rbind(out, data.frame(intervals = interval.end, TTOT = TTOT, events = events, rate = events / TTOT, m2ll = 2 * (rate * TTOT - events * log(rate))))
  }

  return(out)
}
