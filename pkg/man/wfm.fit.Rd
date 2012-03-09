\name{wfm.fit}
\alias{wfm.fit}

\title{Fit Wfm model to trancriptome data}

\description{Main function to fit a wavelet-based functional model to the tiling array expression data.}

\usage{wfm.fit(object, filter.overlap=NULL, design=c("time","circadian","group","factorial","custom"), n.levels, factor.levels=NULL, chromosome, strand, minPos, maxPos, design.matrix=NULL, var.eps=c("margLik","mad"), prior=c("normal","improper"), eqsmooth=FALSE, max.it=20, wave.filt="haar", skiplevels=NULL, trace=FALSE, save.obs=c("plot","regions","all"))}

\arguments{
  \item{object}{object of class \code{WaveTilingFeatureSet}}
  \item{filter.overlap}{object of class \code{mapFilterProbe}}
  \item{design}{character indicating the design of tiling array experiment. Currently, the following designs are implemented: \code{time} for a time-course design based on polynomial contrasts; \code{circadian} for circadian rhythm analysis; \code{group} for unordered one-factor designs; \code{factorial} for two-factor designs; \code{custom} for other designs. When using
  \code{design="custom"} a specific \code{design.matrix} needs to be given.}
  \item{n.levels}{number of levels in wavelet decomposition (integer)}
  \item{factor.levels}{factor levels in case of two-factor analysis. Vector of integers with length the number of factors in the experiment, and with elements the number of levels for the respective factors.}
  \item{chromosome}{numeric to indicate the chromosome associated with transcriptome data to fit}
  \item{strand}{character to indicate the strand orientation associated with transcriptome data to fit. Either \code{"forward"} or \code{"reverse"}.}
  \item{minPos}{integer to indicate minimum genomic position}
  \item{maxPos}{integer to indicate maximum genomic position}
  \item{design.matrix}{custom design matrix to use}
  \item{var.eps}{character indicating how to estimate residual variance. Either \code{"margLik"} for marginal maximum likelihood based estimation or \code{"mad"} for estimation based on the MAD (more info see references).}
  \item{prior}{character indicating which prior distribution to put on effect functions. Either \code{"normal"} for a normally distributed prior, or \code{"improper"} for an improper prior (more info see references).}
  \item{eqsmooth}{logical indicating whether to force equal amount of smooth for different effect functions or not}
  \item{max.it}{integer giving the maximum number of iteration for estimation}
  \item{wave.filt}{character indicating which wavelet filter to use. Default is \code{"haar"}.}
  \item{skiplevels}{integer indicating how many wavelet levels to put equal to 0}
  \item{trace}{logical indicating whether to trace estimation}
  \item{save.obs}{character to indicate which output to store in return object. Either \code{"plot"}: all info needed to make the plots or \code{"all"}: store all possible info.}
}

\value{object of class WfmFit}

\references{
[1] Clement L, De Beuf K, Thas O, Vuylsteke M, Irizarry RA and Crainiceanu CM. (2012) Fast wavelet based functional models for transcriptome analysis with tiling arrays. Statistical Applications in Genetics and Molecular Biology 11: Iss. 1, Article 4.

[2] De Beuf K, Andriankaja, M, Thas O, Inze D, Crainiceanu CM and Clement L (2012) Model-based analysis of tiling array expression studies with flexible designs. Technical document.
}

\author{Kristof De Beuf <kristof.debeuf@ugent.be>}

\examples{
\dontrun{
  library(waveTilingData)
  data(leafdev)
  data(leafdevMapAndFilterTAIR9)
  leafdev <- as(leafdev,"WaveTilingFeatureSet")
  ## background-correction and normalization
  leafdevBQ <- bgCorrQn(leafdev,useMapFilter=leafdevMapAndFilterTAIR9)
  leafdevFit <- wfm.fit(leafdevBQ,filter.overlap=leafdevMapAndFilterTAIR9,design="time",n.levels=10,chromosome=1,strand="forward",var.eps="marg",prior="improper",skiplevels=1,save.obs="plot",trace=TRUE)
}
}

\keyword{manip}
