\name{stat1}
\alias{stat1}
\title{Statistics for one trait}
\description{ 
Calculates the narrow-sense heritability and its confidence intervals for one trait.
}
\usage{
stat1(paraDATA,
      alpha=0.05,
      frommethod="REML",
      Satterthwaite=TRUE,
      tex.table=FALSE,
      file=FALSE,
      path="~/Desktop/myproject.qgen")
}
\arguments{
\item{paraDATA}{object of class \code{\link{paraDATA}} created by the function \code{\link{est}} (or \code{\link{emp}}, see below)}
\item{alpha}{two sided error probability}
\item{frommethod}{If the paraDATA--object provides results form different
  variance--partioning methods you can choose which to take otherwise the REML
  solutions taken without warning:
  \itemize{
    \item{\code{"ANOVAuw"}}{ANOVA estimators with unweighted sums of
      squares}
    \item{\code{"ANOVA"}}{normal ANOVA estimators}
    \item{\code{"REML"}}{REML solutions}
    }}
\item{Satterthwaite}{Logical flag: Should the approximated degree of freedoms for sire and dam be used (only possible if they are available)?}
\item{file}{Logical flag: Should the object be written to \cr
  \code{"~/Desktop/myproject.qgen/stat1.rda"?}}
\item{tex.table}{Logical flag: should a LaTeX file be created with a
  table of all estimates}
\item{path}{the path for the output files}
}
\details{
This function calculates the heritability of individual traits. If the
paraDATA provides parameters for several characters or several
environments or both, then  for every character--trait combination
heritability is calculated individualy. If the paraDATA has a \code{unbal}
slot with data about the balancedness of the data set then also the
confidence intervals are calculated. If the paraDATA provides parameters
estimated with diffent variance partitioning methods then it is possible
to choose among them.}
\value{
  an object of the class \code{\link{stat}}.
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
(not yet implemented) \code{\link{stat2e}}, for the same statistics of an index of two environments; \code{statXt}, for several traits within one environment.
}
\examples{
## calculates narrow-sense heritability and its confidence intervals
myparameters <- est(sim(the()))   # takes the default causal parameters of the()
                                  # simulates a data set
                                  # estimates the variance components
stat1(myparameters)               # calculates heritability and
                                  # confidence invals from estimated parameters
}
\keyword{file}
