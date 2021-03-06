% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/qts.R
\name{add_noise}
\alias{add_noise}
\title{QTS Noise Generator}
\usage{
add_noise(qts, n = 1, alpha = 0.1, beta = 0.5)
}
\arguments{
\item{qts}{A \code{\link[tibble]{tibble}} specifying the input QTS.}

\item{n}{An integer specifying how many noisy QTS should be generated.
Defaults to \code{1L}.}

\item{alpha}{A positive scalar specifying the variance of each component of
the log-QTS. Defaults to \code{0.2}.}

\item{beta}{A positive scalar specifying the exponential weight. Defaults to
\code{0.5}.}
}
\value{
A list of \code{n} QTS with added noise as specified by parameters \code{alpha}
and \code{beta}.
}
\description{
This function adds uncorrelated Gaussian noise to the logarithm QTS using an
exponential covariance function.
}
\details{
See \code{\link[roahd]{exp_cov_function}} for details about the roles of
\code{alpha} and \code{beta} in the definition of the covariance operator.
}
\examples{
qts <- tibble::tibble(
  time = seq(0, 1, len = 101),
  w = 0,
  x = cos(2 * pi * time),
  y = time^2,
  z = sin(time)
)
qts <- exp_qts(qts)
qts_n <- add_noise(qts, n = 10)
}
