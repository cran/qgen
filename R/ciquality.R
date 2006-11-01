ciquality <- function(ciTypeL, ciTypeU, theta){
  Uerror <- NULL; Lerror <- NULL; ciTypeLerror <- NULL; ciTypeUerror <- NULL
  if(!is.null(ciTypeL)&!is.null(ciTypeU)){
    ciTypeLout <- rep(NA, dim(ciTypeL)[2])
    ciTypeLin <- rep(NA, dim(ciTypeL)[2])
    ciTypeLna <-  rep(NA, dim(ciTypeL)[2])
    ciTypeUout <- rep(NA, dim(ciTypeU)[2])
    ciTypeUin <- rep(NA, dim(ciTypeU)[2])
    ciTypeUna <-  rep(NA, dim(ciTypeU)[2])
    den.ciTypeU <- as.list(rep(NA, dim(ciTypeU)[2]))
    den.ciTypeL <- as.list(rep(NA, dim(ciTypeL)[2]))
    ciTypeUout <- colSums(ciTypeU<matrix(theta,ncol=dim(theta)[2],nrow=dim(ciTypeU)[1]),na.rm=TRUE)
    ciTypeUin <- colSums(ciTypeU>matrix(theta,ncol=dim(theta)[2],nrow=dim(ciTypeU)[1]),na.rm=TRUE)
    ciTypeUna <- colSums(is.na(ciTypeU),na.rm=TRUE)
    ciTypeLout <- colSums(ciTypeL>matrix(theta,ncol=dim(theta)[2],nrow=dim(ciTypeL)[1]),na.rm=TRUE)
    ciTypeLin <- colSums(ciTypeL<matrix(theta,ncol=dim(theta)[2],nrow=dim(ciTypeL)[1]),na.rm=TRUE)
    ciTypeLna <- colSums(is.na(ciTypeL),na.rm=TRUE)
      for(cc in 1:dim(ciTypeL)[2]){
          if(ciTypeUout[cc]==dim(ciTypeU)[1]){ciTypeUout[cc] <- NA}
        if(sum(!is.na(ciTypeU[1:dim(ciTypeU)[1],cc]))>1){ #otherwise density is not possible
          den.ciTypeU[[cc]] <- density(ciTypeU[,cc], na.rm=TRUE)
        }
        if(ciTypeLout[cc]==dim(ciTypeL)[1]){ciTypeLout[cc] <- NA}
        if(sum(!is.na(ciTypeL[1:dim(ciTypeL)[1],1]))>1){
          den.ciTypeL[[cc]] <- density(ciTypeL[,cc], na.rm=TRUE)
        }
      }
    ciTypeLerror <- ciTypeLout / (ciTypeLin+ciTypeLout+ciTypeLna)
    ciTypeUerror <- ciTypeUout / (ciTypeUin+ciTypeUout+ciTypeUna)
    ciTypeerror <- ciTypeLerror+ciTypeUerror
    list(Uerror=ciTypeUerror, Lerror=ciTypeLerror, error=ciTypeerror, Udensity=den.ciTypeU, Ldensity=den.ciTypeL)
  }else{
    NULL
  }
}




#ciTypeUout[cc] <- sum(ciTypeU[,cc]<theta[cc])#ciTypeUout[cc] <- length(ciTypeU[ciTypeU[,cc]<theta[cc],cc])
