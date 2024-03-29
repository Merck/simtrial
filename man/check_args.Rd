% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/check_args.R
\name{check_args}
\alias{check_args}
\title{Check argument types, length, or dimension}
\usage{
check_args(arg, type, length = NULL, dim = NULL)
}
\arguments{
\item{arg}{An argument to be checked.}

\item{type}{A character vector of candidate argument type.}

\item{length}{A numeric value of argument length or \code{NULL}.}

\item{dim}{A numeric vector of argument dimension or \code{NULL}.}
}
\value{
Check failure detailed error message.
}
\description{
Check argument types, length, or dimension
}
\details{
If \code{type}, \code{length} or \code{dim} is \code{NULL}, the corresponding check will not be executed.
}
\section{Specification}{

\if{latex}{
 \itemize{
   \item Check if arg is NULL.
   \item Extract the type, length and dim information from arg.
   \item Compare with target values and report error message if it does not match.
 }
 }
\if{html}{The contents of this section are shown in PDF user manual only.}
}

\examples{
\dontrun{
tbl <- as.data.frame(matrix(1:9, nrow = 3))
simtrial:::check_args(arg = tbl, type = c("data.frame"))

vec <- c("a", "b", "c")
simtrial:::check_args(arg = vec, type = c("character"), length = 3)
}
}
\keyword{internal}
