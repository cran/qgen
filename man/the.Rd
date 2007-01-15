\name{the}
\alias{the}
\title{Transformation of causal genetic parameters to observational parameters}
\description{
  For theoretical investigations we need a function that takes the
  chosen causal genetic parameters and transforms them into observational
  parameters. Additionally the form of the object needs to be adjusted
  to fit the form of an object of class \code{\link{paraDATA}}. This function
  handles the following causal variance components: additive; dominance,
  maternal, micro-environmental; and the following sources of observable
  variance components: sire, dam, and individuals.
}
\usage{
the(partitioning = "all",
                chN = 1, # number of characters
                enN = 1, # number of environments
                fbN = 1, # number of fixed blocks
                rbN = 1, # number of random blocks
                siN = 100, # number sires (within each block)
                daN = 6, # number of dams (within each sire)
                idN = 3, # number of individuals (within each dam)
                randomblockCor = matrix(0.5, chN*enN, chN*enN),
                randomblockVar = rep(100, chN*enN),
                additiveCor = matrix(0.5, chN*enN, chN*enN),
                additiveVar = rep(100, chN*enN),
                dominanceCor = matrix(0.5, chN*enN, chN*enN),
                dominanceVar = rep(100, chN*enN),
                maternalCor = matrix(0.5, chN*enN, chN*enN),
                maternalVar = rep(100, chN*enN),
                environmentalCor = matrix(0, chN*enN, chN*enN),
                environmentalVar = rep(100, chN*enN),
                ch.names = paste("ch",1:chN,sep=""),
                en.names = paste("en",1:enN,sep=""),
                fb.names = paste("fb",1:fbN,sep=""),
                fixe = array(0, dim=c(enN, chN, fbN)),
                miss = array(0, dim=c(enN, chN, fbN)),
                file=TRUE,
                path="~/Desktop/myproject.qgen/")
}
\arguments{
  \item{partitioning}{character string: The method of variance
    partitioning used throughout the analysis \describe{
      \item{\code{"all"}}{(default) resampling is based on REML solutions but
	both ANOVA estimators are provided additionally}
      \item{\code{"REML"}}{uses the package \code{lme4} to calculate
	REML solutions}
      \item{\code{"ANOVA"}}{calculates traditional ANOVA estimators; can only
	be used for investigating exactly one character in one environment}
      \item{\code{"ANOVAuw"}}{calculates ANOVA estimators using unweighted sums
	of squares; can only be used for investigating exactly one character
	in one environment}
  }}
  \item{chN}{number of characters}
  \item{enN}{number of environments}
  \item{fbN}{number of fixed blocks}
  \item{rbN}{number of random blocks, a random factor crossed with sire}
  \item{siN}{number of sires}
  \item{daN}{number of dams within sires}
  \item{idN}{number of individuals within dams}
  correlation matrix (environment x character);
    characters are grouped within the environments:
      \item{randomblockCor}{for the random block effects}
      \item{additiveCor}{for the additive genetic effects}
      \item{dominanceCor}{for the dominance genetic effects}
      \item{maternalCor}{for the maternal effects}
      \item{environmentalCor}{for the microenvironmental effects}
  vector with the variance of all environment--character
    combinations; characters are grouped within the environments.
      \item{randomblockVar}{for the random blocks}
      \item{additiveVar}{for the additive genetic effects}
      \item{dominanceVar}{for the dominance genetic effects}
      \item{maternalVar}{for the maternal effects}
      \item{environmentalVar}{for the microenvironmental effects}
  \item{fixe}{array with the mean of every
    fixedblock--environment--character combinations (sequence:
    character in environments in fixedblock, ex. fb1: en1ch1,
    en1ch2...en2ch1...enxchx fb2: en1ch1,...)}
  \item{miss}{vector with the proportion of randomly missing values per environment--character combination; if scalar all are
    assumed to be equal}
  \item{en.names}{vector of character strings: names for the environments}
  \item{ch.names}{vector of character strings: names for the traits}
  \item{fb.names}{vector of character strings: names for the fixed
    blocks}
  \item{file}{logical flag: Should the object be written to a file
    called "the.rda"?}
  \item{path}{character vector: containing a single path name}
}
\details{
This function transforms the causal sources of variance (additive
genetic effects, dominance genetic effects, common (maternal)
environmental effects, and individual (special) environmental effects)
into the three observable variance components from a nested full--sib,
half--sib mating design (North-Carolina Design I). Assumptions: no sources of epistatic variance; see Lynch & Walsh 1997, p.572.
}
\value{
  a \code{\link{paraDATA}} object: with the individual slots
  \item{orig}{fully specified}
  \item{supl}{fully specified}
  \item{para}{fully specified}
  \item{data}{empty}
  \item{spec}{empty}
}
\references{
Lynch, M. and Walsh, B. (1997) \emph{Genetics and analysis of
  quantitative traits.} Sinauer.
}
\seealso{
\code{\link{emp}}, \code{\link{paraDATA}}
}
\note{
\emph{name:} theoretical investigations on chosen biological parameters
}
\examples{
## create a paraDATA-object
## (with a full para section and an empty DATA section)
the()
}
\keyword{file}