\name{sinapis}
\alias{sinapis}
\docType{data}
\title{Paternal Half-sib Data from Sinapis arvensis}
\description{
   Data from an experiment with paternal half-sibs. It allows to study
   the inheritance of leaf area and plant weight and the genetic
   correlation among this two characters.
}
\usage{data(sinapis)}
\format{
  A data frame with 1100 observations on the following 7 variables.
  \describe{
    \item{\code{fixedblock}}{a blocking factor with 3 levels}
    \item{\code{randomblock}}{a blocking factor without any levels}
    \item{\code{sire}}{a factor with 31 levels}
    \item{\code{dam}}{a factor with 3 levels (to be nested within sire)}
    \item{\code{individual}}{a factor with 1100 levels}
    \item{\code{environment}}{a factor with 2 levels ("cont" for
      control, "herb" for herbivory)} 
    \item{\code{leafarea}}{numeric: the area of one leaf in (unit:
      square centimetre)}
    \item{\code{plantweight}}{numeric: the fresh weight of the
      above-ground plant parts (unit: gram)}
  }
}
\source{
 my phd thesis
}
\examples{
data(sinapis)
## str(sinapis)
}
\keyword{datasets}
