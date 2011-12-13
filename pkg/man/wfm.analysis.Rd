\name{wfm.analysis}
\alias{wfm.analysis}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Perform Wfm analysis
}

\description{
## %%  ~~ A concise (1-5 lines) description of what the function does. ~~
}

\usage{
wfm.analysis(object, filter.overlap=NULL, method=c("twoGroup","compareGroupsTime","compareGroupsFactor","circadian","meansByGroupTime","meansByGroupFactor","effectsTime","twoFactors"), n.levels, chromosome, strand, minPos, maxPos, design.matrix=NULL, var.eps=c("margLik","mad"), prior=c("normal","improper"), max.it=20, wave.filt="haar", skiplevels=NULL, trace=FALSE, save.obs=c("plot","regions","all"), alpha=0.05, nsim=1000, delta=NULL, rescale=NULL, two.sided=NULL, minRunPos=90, minRunProbe=1, factor.levels=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
object of class WaveTilingFeatureSet
}
  \item{filter.overlap}{
mapFilterProbe object
}
  \item{method}{
Method to use in Wfm analysis
}
  \item{n.levels}{
number of levels in wavelet decomposition
}
  \item{chromosome}{
chromosome for Wfm analysis
}
  \item{strand}{
strand for Wfm analysis
}
  \item{minPos}{
minimum position
}
  \item{maxPos}{
maximum position
}
  \item{design.matrix}{
custom design matrix to use
}
  \item{var.eps}{
how to estimate residual variance
}
  \item{prior}{
prior distribution on effect functions
}
  \item{max.it}{
maximum number of iteration for estimation
}
  \item{wave.filt}{
wavelet filter to use
}
  \item{skiplevels}{
wavelet levels to put at 0
}
  \item{trace}{
trace estimation
}
  \item{save.obs}{
which output to store in return object
}
  \item{alpha}{
level of significance
}
  \item{nsim}{
number of simulations for inference
}
  \item{delta}{
threshold value
}
  \item{rescale}{
rescale matrix
}
  \item{two.sided}{
one- or two-sided inference
}
  \item{minRunPos}{
minrun by position
}
  \item{minRunProbe}{
minrun by number of probes
}
  \item{factor.levels}{
factor levels in case of two-factor analysis
}
}

\details{
## %%  ~~ If necessary, more details than the description above ~~
}

\value{
object of class Wfm
}

\references{
## %% ~put references to the literature/web site here ~
}

\author{
Kristof De Beuf <kristof.debeuf@ugent.be>
}

\note{
## %%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
## %% ~~objects to See Also as \code{\link{help}}, ~~~
}

\examples{
##
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{manip}
