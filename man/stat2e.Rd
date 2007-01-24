\name{stat2e}
\alias{stat2e}
\title{Statistics for characters in two environments}
\description{ 
Calculates the narrow-sense heritability of plasticity and its confidence intervals.
}
\usage{
stat2e(paraDATA,
       alpha=0.05,
       v=list( c(-1,1),c(0,1),c(1,0),c(1,1),c(1,-1) ),
        Satterthwaite=TRUE,
       file=FALSE,
       path="~/qgen/")
}
\arguments{
\item{paraDATA}{object of class \code{\link{paraDATA}} created by the function \code{\link{est}} (or \code{\link{emp}}, see below)}
\item{alpha}{two sided error probability}
\item{v}{list with weights for the two environments; the indices}
\item{Satterthwaite}{Logical flag: Should the approximated degree of freedoms for sire and dam be used (only possible if they are available)?}
\item{file}{Logical flag: Should the object be written to \cr
  a file called "\code{stat2e.rda}"?}
\item{path}{the path for the output files}
}
\details{
This function calculates the heritability of individual traits. If the
feeder provides parameters for several characters or several
environments or both, then  for every character--trait combination
heritability is calculated individualy. If the feeder has a \code{unbal}
slot with data about the balancedness of the data set then also the
confidence intervals are calculated. If the feeder provides parameters
estimated with diffent variance partitioning methods then it is possible
to choose among them.
}
\value{
  an object of the class \code{\link{stat}}
}
\note{
\emph{Attributes}: \code{limes}
}
\references{
  Burdick, R. K. and Graybill, F. A. (1992)
  \emph{Confidence intervals on variance components.} Marcel Dekker, New
  York.

  Sen, B., Graybill, F. A. and Ting, N. (1992)
  Confidence intervals on variance components. \emph{Biometrical Journal}
  \bold{3}, 259--274.
}
\seealso{
\code{\link{stat1}}, for the same statistics of one trait; \code{statXt}, for several traits within one environment.
}
\examples{
## calculates narrow-sense heritability and its confidence intervals
}
\keyword{file}
