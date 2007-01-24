\name{sta}
\alias{sta}
\title{Calculates the specified statistic on all replications of one
  resampling level}
\description{
This function takes the name of one file and applies the selected
function to all replications of the same resampling level.
}
\usage{
sta(filename= "",
    statistic.name,
    alpha = 0.05,
    transformation = "none",
    file=TRUE,
    path = "~/qgen/")
}
\arguments{
  \item{filename}{character string indicating the name of the file (all
    files of this resampling level will then be chosen!)}
  \item{alpha}{number indicating the two sided error probability}
  \item{statistic.name}{indicate a function that calculates the
    statistic you are interested in (e.g. \code{\link{stat1}})}
  \item{transformation}{character string indicating the transformation
    that should be applied on the statistic to calculate bootstrap
    confidence intervals (\code{"logodds", "asin", "fisher"})}
  \item{file}{logical flag: Should the object be written to a file}
  \item{path}{indicating the path where the files should be read and the
    written}
}
\details{
This function is used to calculate a selected statistic from all
\code{\link{paraDATA}}-objects within the indicated file. It
automatically applies the same statistic to all
\code{\link{paraDATA}}-objects on the same replication level, even if
they are stored in different files.
}
\value{
  \item{multi}{an object of class \code{\link{multi}}.}
  stored in a file within the selected path called "statT.rda",
  "statS.rda", "statR.rda", "statQ.rda", depending on the resampling level.
}
\examples{
##
}
\keyword{file}
