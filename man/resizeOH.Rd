% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/resizeOH.R
\name{resizeOH}
\alias{resizeOH}
\title{Resize onehotTX}
\usage{
resizeOH(onehotTX, window = 250, unify_strand = T)
}
\arguments{
\item{onehotTX}{A \code{GRanges} object whose metadata columns are generated onehotTX.}

\item{window}{An integer indicates the target size which equals to 2*window + 1.}

\item{unify_strand}{A \code{logical} object indicating whether to unify the direction from 5'end to 3'end.}
}
\value{
A \code{data.frame} object contains the resized onehotTX.
}
\description{
This function resizes onehotTX centered on the target site.
}
\examples{
\donttest{
## Resizing a large number of onehotTX is time-consuming.
## Therefore, mcapply from package parallel is suggested.

lapply(onehotTX, resizeOH, window=250)
mcapply(onehotTX, resizeOH, window=250, mc.cores=8)
}
}
