% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ex4_belly.R
\docType{data}
\name{ex4_belly}
\alias{ex4_belly}
\title{Time-to-event data example 4 for non-proportional hazards working group}
\format{
Data frame with 4 variables:
\itemize{
\item \code{id}: Sequential numbering of unique identifiers.
\item \code{month}: Time-to-event.
\item \code{event}: 1 for event, 0 for censored.
\item \code{trt}: 1 for experimental, 0 for control.
}
}
\usage{
data(ex4_belly)
}
\description{
Survival objects reverse-engineered datasets from published Kaplan-Meier
curves.
Individual trials are de-identified since the data are only
approximations of the actual data.
Data are intended to evaluate methods and designs for trials where
non-proportional hazards may be anticipated for outcome data.
}
\examples{
library(survival)

data(ex4_belly)
km1 <- with(ex4_belly, survfit(Surv(month, evntd) ~ trt))
km1
plot(km1)
}
\references{
Lin, Ray S., Ji Lin, Satrajit Roychoudhury, Keaven M. Anderson, Tianle Hu,
Bo Huang, Larry F Leon, Jason J.Z. Liao, Rong Liu, Xiaodong Luo,
Pralay Mukhopadhyay, Rui Qin, Kay Tatsuoka, Xuejing Wang,
Yang Wang, Jian Zhu, Tai-Tsang Chen, Renee Iacona &
Cross-Pharma Non-proportional Hazards Working Group. 2020.
Alternative analysis methods for time to event endpoints under
nonproportional hazards: A comparative analysis.
\emph{Statistics in Biopharmaceutical Research} 12(2): 187--198.
}
\seealso{
\link{ex1_delayed_effect},
\link{ex2_delayed_effect},
\link{ex3_cure_with_ph},
\link{ex5_widening},
\link{ex6_crossing}
}
\keyword{datasets}
