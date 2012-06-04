\name{getSigGenes}
\alias{getSigGenes}

\title{Get significant genes}

\description{Extract the annotated regions (often genes) that overlap with the significant regions found in the wavelet-based transcriptome analysis.}

\usage{getSigGenes(fit,inf,biomartObj)}

\arguments{
  \item{fit}{object of class \code{WfmFit}}
  \item{inf}{object of class \code{WfmInf}}
  \item{biomartObj}{object of class TranscriptDb representing an annotation database generated from BioMart.}
}

\value{GRangesList object. In the elementMetadata of the GRanges elements \code{percOverGene} gives the percentage of basepair overlap of the annotated regions by the detected significant region in the analysis; \code{percOverReg} gives the percentage of basepair overlap of the detected singificant region in the analysis with the annotated region; \code{totPercOverGene} gives per annotated region the total percentage of basepair overlap by all detected significant regions in the analysis that map to that annotated region.}

\author{Kristof De Beuf <kristof.debeuf@ugent.be>}

\examples{
  library(waveTilingData)
  library(TxDb.Athaliana.BioMart.plantsmart12)
  data(leafdevFit)
  data(leafdevInfCompare)
  sigGenesCompare <- getSigGenes(fit=leafdevFit,inf=leafdevInfCompare,	biomartObj=TxDb.Athaliana.BioMart.plantsmart12)
  head(sigGenesCompare[[2]])
}

\keyword{manip}
