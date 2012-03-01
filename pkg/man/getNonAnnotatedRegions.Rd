\name{getNonAnnotatedRegions}
\alias{getNonAnnotatedRegions}

\title{Get non-annotated regions}

\description{Extract the significant regions found in the wavelet-based transcriptome analysis that don't show any overlap with the existing annotation.}

\usage{getNonAnnotatedRegions(fit,inf,annoFile)}

\arguments{
  \item{fit}{object of class \code{WfmFit}}
  \item{inf}{object of class \code{WfmInf}}
  \item{annoFile}{annotation file. This annotation file should be a data.frame with (at least) the following information: \code{chromosome} giving the chromosme of the annotated region; \code{strand} providing the strand orientation of the annotated region; \code{start} and \code{end} containing the start and end position of the annotated region; \code{feature} giving the type of the annotated region (gene, pseudogene, exon,...); \code{ID} containing the annotation ID of the annotated region.}
}

\value{GRangesList object with the non-annotated regions. The first element gives the regions with no annotation overlap on the strand used in the analysis, the second element gives the regions with no annotation overlap on both strands.}

\author{Kristof De Beuf <kristof.debeuf@ugent.be>}

\examples{
\dontrun{
  library(waveTilingData)
  data(tair9gff)
  data(leafdevFit)
  data(leafdevInfCompare)
  nonAnnoCompare <- getNonAnnotatedRegions(fit=leafdevFit,inf=leafdevInfCompare,annoFile=tair9gff)
  nonAnnoCompare
}
}

\keyword{manip}
