\name{filterOverlap}
\alias{filterOverlap}

\title{function to filter probe sequence overlaps and remap probes}

\description{This function remaps the probe sequence to a reference sequence and filters probe sequence overlaps between PM and MM probes and/or between probes on the forward and reverse strand}

\usage{filterOverlap(object, remap = TRUE, BSgenomeObject, chrId, strand = c("forward","reverse","both"), MM=FALSE, ...)}

\arguments{
  \item{object}{object of class \code{WaveTilingFeatureSet}}
  \item{remap}{logical to determine whether the tiling array probe sequences have to be remapped to a more recent reference DNA sequence}
  \item{BSgenomeObject}{object of class \code{BSgenome} containing the genome sequence of the species for which the probes need to be filtered and remapped}
  \item{chrId}{vector of numerics identifying the chromosomes for which the probes have to be filtered and/or remapped}
  \item{strand}{character indicating the strands for which the probes have to be filtered and/or remapped (\code{forward}, \code{reverse} or \code{both})}
  \item{MM}{logical to indicate whether the tiling array contains MM probes or not}
  \item{\dots}{other arguments}
}

\value{An object of class \code{mapFilterProbe} is returned containing the indices of the filtered probes.}

\author{Kristof De Beuf <kristof.debeuf@ugent.be>}

\examples{
\dontrun{
  data(leafdev)
  as(leafdev,"WaveTilingFeatureSet")
  library(BSgenome.Athaliana.TAIR.TAIR9)
  leafdevMapAndFilterTAIR9 <- filterOverlap(leafdev,remap=TRUE,BSgenomeObject=Athaliana,chrId=1:7,strand="both",MM=FALSE)
}
}

\keyword{manip}
