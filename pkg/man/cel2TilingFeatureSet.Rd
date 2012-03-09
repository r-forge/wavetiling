\name{cel2TilingFeatureSet}
\alias{cel2TilingFeatureSet}

\title{Read CEL-files to TilingFeatureSet}

\description{Wrapper function to read in CEL-files and output their content as a \code{TilingFeatureSet}-class object.}

\usage{cel2TilingFeatureSet(dataPath, annotationPackage)}

\arguments{
  \item{dataPath}{character indicating the data path where the CEL-files to read in are stored}
  \item{annotationPackage}{name of the package containing the array probe and annotation information, as produces by the \code{pdInfoBuilder}-package}
}

\details{Uses the \code{list.celfiles} and \code{read.celfiles} functions of the \code{oligo}-package.}

\value{object of class \code{TilingFeatureSet}}

\author{Kristof De Beuf <kristof.debeuf@ugent.be>}

\examples{
  ## No example
}

\keyword{manip}
