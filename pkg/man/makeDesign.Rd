\name{makeDesign}
\alias{makeDesign}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Construct design matrix
}

\description{
## %%  ~~ A concise (1-5 lines) description of what the function does. ~~
}

\usage{
makeDesign(object, method=c("twoGroup","compareGroupsTime","compareGroupsFactor","circadian","twoFactors","meansByGroupTime","meansByGroupFactor","effectsTime"), factor.levels=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
object of class waveTilingFeatureSet
}
  \item{method}{
Wfm method to use
}
  \item{factor.levels}{
Factor levels to use if applying two-factor design
}
}

\details{
## %%  ~~ If necessary, more details than the description above ~~
}

\value{
matrix
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
##%% ~~objects to See Also as \code{\link{help}}, ~~~
}

\examples{
##
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{manip}
