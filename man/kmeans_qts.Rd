% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/kmeans.R
\name{kmeans_qts}
\alias{kmeans_qts}
\title{QTS K-Mean Alignment Algorithm}
\usage{
kmeans_qts(
  qts_list,
  k = 1,
  centroid = "mean",
  diss = "l2",
  warping = "affine",
  iter_max = 20,
  nstart = 1000,
  ncores = 1L
)
}
\arguments{
\item{qts_list}{A list of QTS stored as \code{\link[tibble]{tibble}}s with
columns \code{time}, \code{w}, \code{x}, \code{y} and \code{z}.}

\item{k}{An integer specifying the number of clusters to be formed. Defaults
to \code{1L}.}

\item{centroid}{A string specifying which type of centroid should be used.
Choices are \code{mean} and \code{medoid}. Defaults to \code{mean}.}

\item{diss}{A string specifying which type of dissimilarity should be used.
Choices are \code{l1} and \code{pearson}. Defaults to \code{l2}.}

\item{warping}{A string specifying which class of warping functions should be
used. Choices are \code{none}, \code{shift}, \code{dilation} and \code{affine}. Defaults to
\code{affine}.}

\item{iter_max}{An integer specifying the maximum number of iterations for
terminating the k-mean algorithm. Defautls to \code{20L}.}

\item{nstart}{An integer specifying the number of random restart for making
the k-mean results more robust. Defaults to \code{1000L}.}

\item{ncores}{An integer specifying the number of cores to run the multiple
restarts of the k-mean algorithm in parallel. Defaults to \code{1L}.}
}
\value{
A \code{\link[fdakmapp]{kma}} object storing the results of the best
k-mean alignment algorithm run.
}
\description{
This function massages the input quaternion time series to feed them into the
k-mean alignment algorithm for jointly clustering and aligning the input QTS.
}
\examples{
# TO DO
}
