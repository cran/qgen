\name{cal}
\alias{cal}
\title{
  Organising the resampling procedures.
  }
\description{
  A function that generates a text file containing plain \R code. This
  code  has to be read and run by \R
  (e.g. \code{source()}). Thus a \code{\link{multi}}-object is
  generated from which the statistics can be calculated using the \code{\link{sta}}-function.
}
\usage{
cal(filename, repetitions=2, spliting=1,
    file=TRUE, path="~/qgen/")
}
\arguments{
  \item{filename}{the name of a file containing an object of class
    \code{\link{paraDATA}} or \code{\link{multi}}}
  \item{repetitions}{number of repetitions}
  \item{spliting}{not yet implemented; for spliting tasks among
    different computers}
  \item{file}{logical flag: Should the object be written to a file}
  \item{path}{the path for the files: \code{calfile.r}, and \code{multi...rda}}
}
\details{
The organisation of the resampling procedure with a separate file, which has to be run
line by line, makes the organisation of the resampling very flexible. It allwos
to split up very computer intensive calculations among different
machines. Therefore the file created by this function can simply be
devided into different parts that are then run on different
machines. The simple structure of the \code{\link{multi-class}}--object allows
then easily to put the parts together again.
Another advantage is that errors occuring during optimization are
skipped and don't cause problems.
}
\value{
  a \code{\link{multi}} object, stored in a file: path/\code{multi...rda}
}
\references{
}
\seealso{
The vignette and the functions \code{\link{sta}} and  \code{\link{dis}}
}
\note{
\emph{name:} from caling the different functions
}
\keyword{file} 