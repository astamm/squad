% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{resample_qts}
\alias{resample_qts}
\title{QTS Resampling}
\usage{
resample_qts(
  qts,
  tmin = NA_real_,
  tmax = NA_real_,
  nout = 0L,
  disable_normalization = FALSE
)
}
\arguments{
\item{qts}{A quaternion time series stored as a \code{\link[tibble]{tibble}}
with columns \code{time}, \code{w}, \code{x}, \code{y} and \code{z}.}

\item{tmin}{A numeric value specifying the lower bound of the time interval
over which uniform resampling should take place. It must satisfy \code{tmin >= min(qts$time)}. Defaults to \code{NA} in which case it is set to
\code{min(qts$time)}.}

\item{tmax}{A numeric value specifying the upper bound of the time interval
over which uniform resampling should take place. It must satisfy \code{tmax <= max(qts$time)}. Defaults to \code{NA} in which case it is set to
\code{max(qts$time)}.}

\item{nout}{An integer specifying the size of the uniform grid for time
resampling. Defaults to \code{0L} in which case it uses the same grid size as
the input QTS.}

\item{disable_normalization}{A boolean specifying whether quaternion
normalization should be disabled. Defaults to \code{FALSE} in which case the
function makes sure that quaternions are normalized prior to performing
SLERP interpolation.}
}
\value{
A quaternion time series stored as a \code{\link[tibble]{tibble}}
with columns \code{time}, \code{w}, \code{x}, \code{y} and \code{z} in which quaternions are
uniformly sampled in the range \verb{[tmin, tmax]}.
}
\description{
This function performs uniform resampling using SLERP.
}
\examples{
# TO DO
}
