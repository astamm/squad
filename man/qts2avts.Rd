% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{qts2avts}
\alias{qts2avts}
\title{QTS Transformation to Angular Velocity Time Series}
\usage{
qts2avts(qts, fixed_frame = "global")
}
\arguments{
\item{qts}{A QTS stored as a \code{\link[tibble]{tibble}}s with columns
\code{time}, \code{w}, \code{x}, \code{y} and \code{z}.}

\item{fixed_frame}{A string specifying the fixed frame with respect to which
coordinates of the angular velocity should be computed. Choices are
\code{"global"} or \code{"body"}. Defaults to \code{"global"}.}
}
\value{
A time series stored as a \code{\link[tibble]{tibble}} with columns
\code{time}, \code{x}, \code{y} and \code{z} containing the angular velocity at each time
point.
}
\description{
This function projects a quaternion time series into the space of angular
velocities.
}
\examples{
# TO DO
}
