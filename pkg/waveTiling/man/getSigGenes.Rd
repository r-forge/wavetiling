\name{getSigGenes}
\alias{getSigGenes}

\title{Get significant genes}

\description{Extract the annotated regions (often genes) that overlap with the significant regions found in the wavelet-based transcriptome analysis.}

\usage{getSigGenes(fit,inf,annoFile)}

\arguments{
  \item{fit}{object of class \code{WfmFit}}
  \item{inf}{object of class \code{WfmInf}}
  \item{annoFile}{annotation file. This annotation file should be a data.frame with (at least) the following information: \code{chromosome} giving the chromosme of the annotated region; \code{strand} providing the strand orientation of the annotated region; \code{start} and \code{end} containing the start and end position of the annotated region; \code{feature} giving the type of the annotated region (gene, pseudogene, exon,...); \code{ID} containing the annotation ID of the annotated region.}
}

\value{GRangesList object. In the elementMetadata of the GRanges elements \code{percOverGene} gives the percentage of basepair overlap of the annotated regions by the detected significant region in the analysis; \code{percOverReg} gives the percentage of basepair overlap of the detected singificant region in the analysis with the annotated region; \code{totPercOverGene} gives per annotated region the total percentage of basepair overlap by all detected significant regions in the analysis that map to that annotated region.}

\author{Kristof De Beuf <kristof.debeuf@ugent.be>}

\examples{
\dontrun{
  library(waveTilingData)
  data(tair9gff)
  data(leafdevFit)
  data(leafdevInfCompare)
  sigGenesCompare <- getSigGenes(fit=leafdevFit,inf=leafdevInfCompare,annoFile=tair9gff)
  head(sigGenesCompare[[2]])
}
}

\keyword{manip}
