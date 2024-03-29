% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/scUtilities.R
\name{GetTSSPositions}
\alias{GetTSSPositions}
\title{Find transcriptional start sites}
\usage{
GetTSSPositions(ranges, biotypes = "protein_coding")
}
\arguments{
\item{ranges}{A GRanges object containing gene annotations.}

\item{biotypes}{Gene biotypes to include. If NULL, use all biotypes in the
supplied gene annotation.}
}
\value{
transcriptional start sites
}
\description{
Find transcriptional start sites
}
\details{
Get the TSS positions from a set of genomic ranges containing gene positions.
Ranges can contain exons, introns, UTRs, etc, rather than the whole
transcript. Only protein coding gene biotypes are included in output.
}
\keyword{internal}
