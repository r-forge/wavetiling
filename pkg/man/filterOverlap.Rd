\name{filterOverlap}
\alias{filterOverlap}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
function to filter probe sequence overlaps
}
\description{
This function remaps the probe sequence to a reference sequence and filters probe sequence overlaps between PM and MM probes and/or between probes on the forward and reverse strand
}
\usage{
filterOverlap(object, remap = TRUE, fastaFile, chrId, strand = c("forward","reverse","both"), MM=FALSE, ...)
}

\arguments{
  \item{object}{
object of class waveTilingFeatureSet
}
  \item{remap}{
determines whether the tiling array probe sequences have to be remapped to a more recent reference DNA sequence
}
  \item{fastaFile}{
fasta file containing the reference DNA sequence to remap the tiling array probes against
}
  \item{chrId}{
identifies the chromosome for the tiling array under consideration
}
  \item{strand}{
identifies the corresponding strand for the tiling array under consideration
}
  \item{MM}{
does the tiling array contain MM probes?
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
}
}
\details{
## %%  ~~ If necessary, more details than the description above ~~
}

\value{
Returns an object of class filterOverlap
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
## to do
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{manip}
