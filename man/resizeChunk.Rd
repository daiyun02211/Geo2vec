% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/resizeChunk.R
\name{resizeChunk}
\alias{resizeChunk}
\title{Resize chunkTX.}
\usage{
resizeChunk(chunkTX, window = 17, unify_strand = TRUE, region = FALSE)
}
\arguments{
\item{chunkTX}{A \code{GRanges} object whose metadata columns are generated chunkTX.}

\item{window}{An integer indicates the target size which equals to 2*window + 1.}

\item{unify_strand}{A \code{logical} object indicating whether to unify the direction from 5'end to 3'end.}

\item{region}{A \code{logical} object indicating whether to resize region chunkTX.}
}
\value{
A \code{data.frame} object contains the resized chunkTX.
}
\description{
This function resizes chunkTX centered on the target site/region.
}
\examples{
\donttest{
## Resizing a large number of chunkTX is time-consuming.
## Therefore, mcapply from package parallel is suggested.

lapply(chunkTX, resizeChunk, window=17, region=FALSE)
mcapply(chunkTX, resizeChunk, window=17, region=FALSE, mc.cores=8)
}
}
