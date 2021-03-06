% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/generics.R
\name{summary.ThinData}
\alias{summary.ThinData}
\title{Provide summary output of a ThinData S3 object.}
\usage{
\method{summary}{ThinData}(object, ...)
}
\arguments{
\item{object}{A ThinData S3 object. This is generally output by either
\code{\link{thin_diff}}, \code{\link{thin_2group}}, or
\code{\link{thin_lib}}.}

\item{...}{Not used.}
}
\value{
Returns nothing. Prints out some summary information on
    \code{object}.
}
\description{
Provide summary output of a ThinData S3 object.
}
\author{
David Gerard
}
