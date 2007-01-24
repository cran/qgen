\name{emp}
\alias{emp}
\title{Transforms empirical data sets into a conformable paraDATA object}
\description{ 
Takes empirical data from a dataframe, defines
new columns, and puts the data in an object of the \code{\link{paraDATA-class}}.
}
\usage{
emp(data.use=NULL,
    fixedblock.use, character.use,
    environment.use="environment",
    randomblock.use,
    sire.use="sire",
    dam.use="dam",
    individual.use="individual",
    without="", partitioning="REML",
    file=TRUE, path="~/qgen/")
}
\arguments{
\item{data.use}{the name of a dataframe}
\item{fixedblock.use}{the name of the column in \code{data.use}
containing the levels of a fixed blocking factor (as factor)}
\item{character.use}{a vector of characters strings indicating the column
  names to be chosen for analysis (numeric values only)}
\item{environment.use}{the name of the column in \code{data.use}
containing the levels of the environments (as factor)}
\item{randomblock.use}{the name of the column in \code{data.use}
indicating the levels of a random blocking factor (as factor)}
\item{sire.use}{the name of the column in \code{data.use} containing the
levels of the sire effect (as factor)}
\item{dam.use}{the name of the column in \code{data.use} containing the
levels of the dam effect (as factor)}
\item{individual.use}{the name of the column in \code{data.use} containing the
levels of the individual effect (as factor)}
\item{without}{a vector of characters indicating the names of the rows
  to be excluded (for outliers)}
\item{partitioning}{a character string indicating the variance
partitioning method (used by the function \code{\link{est}}); see
\code{\link{the}} for details}
\item{file}{logical flag: Should the object be written to \code{emp.rda}}
\item{path}{indicating the path where the file \code{emp.rda} should be written}
}
\value{an object of class \code{\link{paraDATA}} with
  \item{\link{orig}}{ \describe{
      \item{hist}{"emp"}
  }}
}
\seealso{
\code{\link{paraDATA}}
}
\note{
  \emph{name:} empirical investigations
}
\examples{
emp(file=FALSE) #creates a paraDATA object with the data from "sinapis"
}
\keyword{file}