% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/milestone.R
\name{milestone}
\alias{milestone}
\title{Milestone test for two survival curves}
\usage{
milestone(data, ms_time)
}
\arguments{
\item{data}{Data frame containing at least 3 columns:
\itemize{
\item \code{tte} - Time to event.
\item \code{event} - Event indicator.
\item \code{treatment} - Grouping variable.
}}

\item{ms_time}{Milestone analysis time.}
}
\value{
A data frame containing:
\itemize{
\item \code{method} - The method, always \code{"milestone"}.
\item \code{z} - Test statistics.
\item \code{ms_time} - Milestone time point.
\item \code{surv_ctrl} - Survival rate of the control arm.
\item \code{surv_exp} - Survival rate of the experimental arm.
\item \code{surv_diff} - Survival difference between the experimental and control arm.
\item \code{std_err_ctrl} - Standard error of the control arm.
\item \code{std_err_exp} - Standard error of the experimental arm.
}
}
\description{
Milestone test for two survival curves
}
\examples{
sim_pw_surv(n = 200) |>
  cut_data_by_event(150) |>
  milestone(10)
}
\references{
Klein, J. P., Logan, B., Harhoff, M., & Andersen, P. K. (2007).
"Analyzing survival curves at a fixed point in time."
\emph{Statistics in Medicine}, 26(24), 4505--4519.
}
