\name{qgen-classes}
\alias{paraDATA-class}
\alias{orig-class}
\alias{supl-class}
\alias{para-class}
\alias{DATA-class}
\alias{spec-class}
\alias{multi-class}
\alias{stat-class}
\alias{paraDATA}
\alias{orig}
\alias{supl}
\alias{para}
\alias{DATA}
\alias{spec}
\alias{multi}
\alias{stat}
\title{qgen: Class definitions}
\description{Describing classes associated with qgen}
\usage{
###
"paraDATA"# is a virtual class, extended by:
"orig"
"supl"
"para"
"DATA"
"spec"
###
"multi"
###
"stat"
}
\section{paraDATA}{\describe{
    \item{orig}{\describe{
	\item{hist}{A list with the name of all functions that have
	  handled the object; in chronological order}
	\item{warn}{A list with warnings (character strings) from the functions that have
	  handled the object; in chronological order}
	\item{time}{A list with the time  (character strings) when the
	  functions have handled the object; in chronological order}
	\item{part}{A character string; defines which variance
	  partitioning method(s) should be used: \describe{
	    \item{\code{"all"}}{REML for the resampling; additionally ANOVA and
	      ANOVAuw}
	    \item{\code{"REML"}}{REML}
	    \item{\code{"ANOVA"}}{ANOVA with traditional sums of squares}
	    \item{\code{"ANOVAuw"}}{ANOVA with unweighted sums of squares}
	}}
    }}
    \item{supl}{\describe{
	\item{chN}{the number of characters}
	\item{enN}{the number of environments}
	\item{fbN}{the number of fixed blocks (crossed with
	  environments and  characters)}
	\item{rbN}{the number of random blocks (crossed with sire)}
	\item{siN}{the number of sires}
	\item{daN}{the number of dams within sires}
	\item{idN}{the number of individuals}
	\item{miss}{proportion of randomly missing values per
	  environment--character--fixedblock combination in a array; \code{dim=c(fb,en,ch)}; if scalar all are
	  assumed to be equal}
    }}
    \item{para}{all parameters needed to describe the distribution of
      the data:\describe{
	\item{rbS}{matrix, random block effect}
	\item{siS}{matrix, sire effect}
	\item{daS}{matrix, dam effect}
	\item{idS}{matrix, individual effect}
	\item{phs}{matrix, phenotypic effect (not observable!)}
	\item{error}{scalar, variance of model residuals}
	\item{fixe}{array  (\code{dim=c(fb,en,ch)}), cell means of every environment--characer-fixedblock
	  combination}
    }}
    \item{DATA}{\describe{
	\item{dat}{data-frame with columns \describe{
	    \item{\code{ch}}{a factor with one level for every character}
	    \item{\code{en}}{a factor with one level for every environment}
	    \item{\code{fb}}{a factor with one level for every fixed block
	      (crossed with the characters and environment)}
	    \item{\code{rb}}{a factor with one level for every random block
	      (crossed with the characters and environment)}
	    \item{\code{si}}{a factor with one level for every sire}
	    \item{\code{da}}{a factor with one level for every dam}
	    \item{\code{id}}{a factor with one level for every individual}
	    \item{\code{y}}{measurement (numeric)}
    }}}}
    \item{spec}{the place for specific information on the model used to
      estimate the parameter etc. \describe{
      \item{additional.partitioning}{A list with qgen \code{para-class}
	objects containing parameters variance partitioning by ANOVA and/or ANOVAuw (with unweighted sums of
	squares)}
      \item{unbalanced}{A list; information on the unbalancedness of the
	data (where "xx" is "si", "da", or "id")\describe{
	  \item{\code{xxDF}}{normal degrees of freedom}
	  \item{\code{xxDFappREML}}{Satterthwaite approximated degrees of
	    freedom using the REML variance components}
	  \item{\code{xxDFappANOVA}}{Satterthwaite approximated degrees of
	    freedom using the ANOVA estimators for variance
	    components}
	  \item{\code{xxDFappANOVAuw}}{Satterthwaite approximated degrees of
	    freedom using the ANOVA estimators with unweighted sums
	    of squares for variance components}
	}
	weights used to calculate the expected mean squares \describe{
	  \item{\code{w1u}}{for dam variance in calculation of sire mean square}
	  \item{\code{w2u}}{for sire variance}
	  \item{\code{w3u}}{for dam variance in calculation of dam mean square}
      }}
      \item{modelsummary}{A list with additional information on the REML model that was
	used to estimate the parameters \describe{
	  \item{\code{sBLUP}}{A matrix with the best linear unbiased predictors for
	    the sire effect}
	  \item{\code{modelsumry}}{An object of
	    \code{summary.lmer-class} (package Matrix)}
      }}
      \item{secondcontrast}{A list with additional modelsummaries if
	different model contrasts were used}
  }}
}}
\section{multi}{\describe{
    \item{list.paraDATA}{a list of \code{\link{paraDATA}}--objects}
    \item{level}{a character indicating the level of resampling
      (T,S,R,Q)}
    \item{x}{still empty}
    \item{y}{still empty}
}}
\section{stat}{\describe{
    \item{orig}{Object of \code{\link{orig-class}}}
    \item{stat}{A numeric value of the statistic}
    \item{lower.ci}{A numeric value, the lower confidence limit}
    \item{upper.ci}{A numeric value, the upper confidence limit}
    \item{lower.limes}{A numeric value, the lowest possible value of that statistic (for
      setting limits in plots)}
    \item{upper.limes}{A numeric value, the highest possible value of that statistic (for
      setting limits in plots)}
}}
\seealso{
  Functions that construct \code{paraDATA}--objects:
  \code{\link{emp}} (from empirical data), \code{\link{the}} (from
  chosen biological parameters). Functions that handle and manipulate \code{paraDATA}--objects:
  \code{\link{sim}}, \code{\link{est}}, \code{\link{stat1}}.

  The function that constructs \code{multi}--objects:
  \code{\link{cal}} (exactly: running the calfile.rda constructed by \code{\link{cal}}).
  
  The function that handles and manipulates \code{multi}--objects:
  \code{\link{sta}}.

  The function that constructs \code{stat}--objects:
  \code{\link{sta}}.
}
\keyword{}