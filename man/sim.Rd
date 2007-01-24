\name{sim}
\alias{sim}
\title{simulation of quantitative genetic data sets}
\description{
  This function simulates data sets according to the sampling design and
  the parameters given by a \code{\link{paraDATA}}-object.
  }
\usage{
sim(paraDATA, file=FALSE, path="~/qgen/")
}
\arguments{
  \item{paraDATA}{an object of class \code{\link{paraDATA}}}
  \item{file}{Logical flag: Should the \code{\link{paraDATA}} object be
    written to path/"sim.rda"?}
  \item{path}{indicating the path where the file \code{sim.rda} should be written}
}
\value{a \code{\link{paraDATA}}--object, with the slots:
  \item{orig}{full}
  \item{supl}{full}
  \item{para}{full, the realized sampling variance-covariance
    components; \code{phS=NULL}}.
  \item{DATA}{full; the simulated data set}
  \item{spec}{NULL}
}
\note{
  \emph{Name:} simulation
}
\examples{
   my.parameters <- the(file=FALSE)         # defines a paraDATA object
   my.paraDATA <- sim(my.parameters)        # simulates the corresponding data set
   data.set <- my.paraDATA@DATA@dat         # extracts the simulated data
                                            # from the paraDATA--object
}
\keyword{file}