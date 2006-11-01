\name{est}
\alias{est}
\title{estimating parameters from data sets}
\description{
  estimates variance components for random effects and coefficients
  of the fixed effects for a given \code{\link{paraDATA}}--object.
}
\usage{est(paraDATA, file=TRUE, path="~/Desktop/myproject.qgen/")
}
\arguments{
\item{paraDATA}{list with a definded structure; see \code{\link{paraDATA}}}
\item{file}{Logical flag: Should the \code{\link{paraDATA}} object be written to
  path/"est.rda"?}
\item{path}{indicating the path where the file \code{est.rda} should be written}
}
\details{The parameters of the \code{\link{the}} and the
  \code{\link{emp}} function determine the estimations that are performed and
  the exact structure of the output}
}
\value{a \code{\link{paraDATA}}--object:
  \item{orig}{full}
  \item{supl}{full}
  \item{para}{full, with \code{phS=NULL}}
  \item{spec}{depending on the \code{supl}-slot of the paraDATA and on the
    \code{hist}ory of the object itself:\describe{
      \item{unbal}{information on the unbalacedness of the
	data; different degrees of freedom: \describe{
	  \item{xxDF}{normal degrees of freedom}
	  \item{xxDFappREML}{Satterthwaite approximated degrees of
	    freedom using the REML variance components}
	  \item{xxDFappANOVA}{Satterthwaite approximated degrees of
	    freedom using the ANOVA estimators for variance
	    components}
	  \item{xxDFappANOVAuw}{Satterthwaite approximated degrees of
	    freedom using the ANOVA estimators with unweighted sums
	    of squares for variance components}
	}
	weights used to calculate the expected mean squares \describe{
	  \item{\code{w1u}}{for dam variance in calculation of sire mean square}
	  \item{\code{w2u}}{for sire variance}
	  \item{\code{w3u}}{for dam variance in calculation of dam mean square}
      }}
      \item{modelsummary}{additional information:\describe{
	  \item{\code{sBLUP}}{the best linear unbiased predictors for the sire
	    effect}
	  \item{\code{modelsumry}}{the summary(lme4-object)}
      }}
      \item{secondcontrast}{the same model with a different contrast
	matrix to estimate the difference in on character between
	environments: \describe{
	  \item{\code{FIX}}{an array with the means of all fixed effects}
	  \item{\code{modelsumry2}}{the summary(lme4-object) for a second
	    model estimated with contrast for testing whether a character
	    differs between two environments}  
      }}
      \item{\code{secondpartitioning}}{variance components
	additionally (to REML) estimated by:\describe{
	  \item{\code{ANOVA}}{ estimators (normal): \describe{
	      \item{\code{rbS}}{variance component for factor sire}
	      \item{\code{siS}}{variance component for factor sire}
	      \item{\code{daS}}{variance component for factor dam}
	      \item{\code{idS}}{variance component for residuals}
	      \item{\code{phS}}{\code{NULL}}
	      \item{\code{error}}{0; the model residuals are included in the idS}
	  }}
	  \item{\code{ANOVAuw}}{ estimators based on unweighted sums of
	    squares (for balanced data equal to ANOVA): \describe{
	      \item{\code{rbS}}{variance component for factor sire}
	      \item{\code{siS}}{variance component for factor sire}
	      \item{\code{daS}}{variance component for factor dam}
	      \item{\code{idS}}{variance component for residuals}
	      \item{\code{phS}}{\code{NULL}}
	      \item{\code{error}}{0; the model residuals are included in the idS}
	  }}
      }}
  }}
}
\note{
  \emph{name:}est from estimation
}
\references{
  Burdick, R. K. and Graybill, F. A. (1992)
  \emph{Confidence intervals on variance components.} Marcel Dekker.

  Sen, B., Graybill, F. A. and Ting, N. (1992)
  Confidence intervals on variance components. \emph{Biometrical Journal}
  \bold{3}, 259--274.
}
\seealso{
\code{\link{paraDATA}}
}
\examples{
parameters <- est(sim(the()))
   # takes the default causal parameters of the()
   # simulates a data set
   # estimates the parameters
}
\keyword{file}