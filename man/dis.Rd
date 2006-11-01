\name{dis}
\alias{dis}
\title{
  Bootstrap confidence intervals
  }
\description{
  Calculates different bootstrap confidence intervals.
}
\usage{
dis(path="~/Desktop/myproject.qgen/", alpha=0.05)
}
\arguments{
  \item{path}{path searched for \code{stat}X\code{.rda}-files.}
  \item{alpha}{number indicating the two sided error probability}
 }
\details{
Depending on the available levels of resampling, percentile , basic (for
simple resampling) and studentized (for nested resampling) confidence
intervals are calculated
}
\value{
  the confidence intervals are printed
}
\examples{
### dis()
}
\keyword{file}