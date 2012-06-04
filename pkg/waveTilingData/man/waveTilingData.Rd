\name{waveTilingData-package}
\alias{waveTilingData-package}
\alias{waveTilingData}
\docType{package}
\title{
Example datasets for the waveTiling package
}
\description{
This package contains the datasets used in the waveTiling package vignette and examples.
}
\seealso{\code{\link[waveTilingData:leafdev]{leafdev}},\code{\link[waveTilingData:leafdevBQ]{leafdevBQ}},\code{\link[waveTilingData:leafdevFit]{leafdevFit}},\code{\link[waveTilingData:leafdevInfCompare]{leafdevInfCompare}},\code{\link[waveTilingData:leafdevMapAndFilterTAIR9]{leafdevMapAndFilterTAIR9}}
}
\author{
Kristof De Beuf <kristof.debeuf@ugent.be>
}
\examples{
dataDir <- system.file("data", package="waveTilingData")
setwd(dataDir)
dir()
}
