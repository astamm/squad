% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{median_qts}
\alias{median_qts}
\title{QTS Geometric Median}
\usage{
median_qts(qts_list)
}
\arguments{
\item{qts_list}{A list of quaternion time series stored as a
\code{\link[tibble]{tibble}}s with columns \code{time}, \code{w}, \code{x}, \code{y} and \code{z}.}
}
\value{
A quaternion time series stored as a \code{\link[tibble]{tibble}}
with columns \code{time}, \code{w}, \code{x}, \code{y} and \code{z} in which quaternions are the
pointwise geometric median.
}
\description{
This function computes the pointwise geometric median of a list of QTS.
}
\examples{
# TO DO
}
