colVars <- function(x, na.rm=FALSE, dims=1, unbiased=TRUE, SumSquares=FALSE, 
  twopass=FALSE) { 
  if (SumSquares) return(colSums(x^2, na.rm, dims)) 
  N <- colSums(!is.na(x), FALSE, dims) 
  Nm1 <- if (unbiased) N-1 else N 
  if (twopass) {x <- if (dims==length(dim(x))) x - mean(x, na.rm=na.rm) else 
  sweep(x, (dims+1):length(dim(x)), colMeans(x,na.rm,dims))} 
  (colSums(x^2, na.rm, dims) - colSums(x, na.rm, dims)^2/N) / Nm1 
} 
# ------------------------------------------------------------------------------- 
#With twopass=TRUE, this should act very much like the S-Plus function by the 
# same name. (twopass=TRUE is needed for numerical accuracy if mu >> sigma.) 

#For a matrix (as opposed to a >2-dimensional array), with unbiased=FALSE, you'd 
# get the same result as Thomas Lumley's suggestion: 
#   R> colMeans(m*m)-colMeans(m)^2 
# However, this is *not* the same as apply(m, 2, var) because (sample) variances 
# are normally calculated with unbiased=TRUE (i.e. with N-1 in the denominator). 
# IF you have no NA's, then I'd modify Thomas's suggestion slightly: 
#   R> N <- nrow(m) 
#   R> N/(N-1) * (colMeans(m*m)-colMeans(m)^2) 

# hab ich aus dem Internet
