\name{qgen-internal}
\alias{colVars}
\alias{interpolation}
\alias{leading}
\alias{ciquality}
\alias{chiplot}
\alias{interpolation}
\alias{matrixscale}
\alias{pcorr}
\title{Internal qgen Functions}
\description{
  Internal Qgen functions:
  ciquality: A function to evaluate the empirical error rate of confidence intervals.
  chiplot: Plots a diagnostic chiplot
  interpolation: linear interpolation if \eqn{(R+1)\alpha}{(R+1)alpha}
  is not an integer
  leading: A function that takes an integer (or a vector of integers) and
  returns a character string (or vector of character strings) with the
  constan length of 4; with leading zeros
  colVars: Form column variance for numeric arrays
  matrixscale: Scales a matrix, which can be a correlation or variance-covariance matrix
  pcorr: Takes a correlation matrix and returns a partial correlation matrix
}
\usage{
ciquality(ciTypeL, ciTypeU, theta)
chiplot(y, id = FALSE)
interpolation(R, alpha, sort.mat)
leading(x)
colVars(x, na.rm=FALSE, dims=1, unbiased=TRUE,
        SumSquares=FALSE, twopass=FALSE)
matrixscale(Sigma)
pcorr(corr)
}
\arguments{
  \item{ciTypeL}{vector of lower confidence limits}
  \item{ciTypeU}{vector of upper confidence limits}
  \item{theta}{the known value for the statisitc}
  \item{y}{numeric vector}
  \item{id}{logical flag: identifying the individual values in the plot}
  \item{R}{scalar, the number of resamples}
  \item{alpha}{a single number or a vector; indicating the two sided
    error probability}
  \item{sort.mat}{an array with the first dimensions beeing repetitions
    at the same level, the second dimensions are the sorted values used
    for interpolation}
  \item{x}{an integer; if a numeric value is given, as.integer() is
    used}
  
  \item{x}{A numeric array (or a dataframe to convert to a matrix).} 
  \item{na.rm}{Logical: Remove NA's?} 
  \item{dims}{Number of dimensions to sum over [colSums] or leave alone 
              [rowSums]. Only useful when x is a multidimensional array.} 
  \item{unbiased}{Logical: Use (N-1) in the denominator when calculating 
                  variance?} 
  \item{SumSquares}{Logical: If TRUE, colVars just returns sums of squares.} 
  \item{twopass}{Logical: If TRUE, colVars uses the corrected two-pass 
     algorithm of Chan Golub & LeVeque, which is slower but less subject 
     to roundoff error.}
   \item{Sigma}{the matrix to be scaled}
   \item{corr}{a correlation matrix}
  } 
\details{
These are not to be called by the user (or in some cases are just
waiting for proper documentation to be written ;-).
colVars() is a suggestion from the R-mailing list.
}
\author{colVars: Originally by Douglas Bates <bates@stat.wisc.edu> as package 
       "MatUtils". Modified, expanded, and renamed by David Brahm 
       <brahm@alum.mit.edu>, with help of course from the R-help gurus.} 
}
\references{
  Davison, A. C. \& Hinkley, D. V. (1997) \emph{Bootstrap methods and their 
    application.} Cambridge University Press.
}
\keyword{internal}
