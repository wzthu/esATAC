% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/scUtilities.R
\name{StringToGRanges}
\alias{StringToGRanges}
\title{StringToGRanges}
\usage{
StringToGRanges(regions, sep = c("-", "-"), ...)
}
\arguments{
\item{regions}{Vector of genomic region strings}

\item{sep}{Vector of separators to use for genomic string. First element is
used to separate chromosome and coordinates, second separator is used to
separate start and end coordinates.}

\item{...}{Additional arguments passed to
\code{\link[GenomicRanges]{makeGRangesFromDataFrame}}}
}
\value{
Returns a GRanges object
}
\description{
String to GRanges
}
\details{
Convert a genomic coordinate string to a GRanges object
}
\keyword{internal}
