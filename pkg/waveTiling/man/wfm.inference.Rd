\name{wfm.inference}
\alias{wfm.inference}

\title{Perform transcriptome analysis on fitted wavelet-based functional model}

\description{Main function to perform trancriptome analysis on a fitted wavelet-based functional model of class \code{WfmFit}.}

\usage{wfm.inference(object, contrast.matrix=NULL, contrasts=c("compare","means","effects","overallMean"), delta=NULL, two.sided=NULL, minRunPos=90, minRunProbe=1, alpha=0.05, nsim=1000, rescale=NULL)}

\arguments{
  \item{object}{object of class \code{WfmFit}}
  \item{contrast.matrix}{custom contrasts matrix}
  \item{contrasts}{character indicating the type of transcriptome analysis is to be applied. Currently the following types are implemented: \code{compare} for doing a pairwise differential expression analysis between any combination of two groups; \code{effects}, which corresponds to a circadian rhythm analysis if a circadian design is used for the fit, and to a time effects analysis (linear, quadratic,...) if a time-course design is used for the fit; \code{means} for doing a group-wise transcript discovery analysis.}
  \item{delta}{threshold value to be used in the inference procedure. This should be a numeric vector with as first element the threshold for the overall mean transcript discovery and the other elements the threshold for the differential expression, the effects analysis or group-wise mean analysis. If the threshold should be equal for all comparisons, effects or group-wise means only a vector of length 2 is needed. Otherwise, the vector must be of length \code{r+1} with \code{r} the number of pairwise comparisons, effects or group-wise means.}
  \item{two.sided}{logical indicating if one-sided or two-sided tests are desired}
  \item{minRunPos}{minrun by position. An integer to indicate the minimum number of basepairs a significant genomic region should contain.}
  \item{minRunProbe}{minrun by probes. An integer to indicate the minimum number of probes the significant genomic region should map to.}
  \item{alpha}{significance level}
  \item{nsim}{number of simulations used when doing circadian rhythm inference}
  \item{rescale}{rescale matrix}
}

\value{object of class \code{WfmFit}}

\references{
[1] Clement L, De Beuf K, Thas O, Vuylsteke M, Irizarry RA and Crainiceanu CM. (2012) Fast wavelet based functional models for transcriptome analysis with tiling arrays. Statistical Applications in Genetics and Molecular Biology 11: Iss. 1, Article 4.

[2] De Beuf K, Andriankaja, M, Thas O, Inze D, Crainiceanu CM and Clement L (2012) Model-based analysis of tiling array expression studies with flexible designs. Technical document.
}

\author{Kristof De Beuf <kristof.debeuf@ugent.be>}

\examples{
  library(waveTilingData)
  data(leafdevFit)
  delta <-  log(1.2,2)
  leafdevInfCompare <- wfm.inference(leafdevFit,contrasts="compare",delta=c("median",delta))
}

\keyword{manip}
