matrixscale <- function(Sigma){
  if (dim(Sigma)[1]==dim(Sigma)[2]){
    diag(Sigma)[diag(Sigma)<0] <- NA
    n <- dim(Sigma)[1]
    denominator <- (matrix(sqrt(diag(Sigma)),nrow=n, ncol=n) * matrix(sqrt(diag(Sigma)), nrow=n, ncol=n, byrow=TRUE))
    denominator[denominator==0] <- NA
    Corr <- Sigma / denominator
    Corr
  }else{matrix(NA, nrow=dim(Sigma)[1], ncol=dim(Sigma)[2])}
}

## Scales a Matrix with c_ij<-s_ij/sqrt(s_ii*s_jj)
## if one variance is negative or excatly equal to zero it indicates the corresponding elements of the correlation matrix with a NA

## a better ERROR message is needed !
