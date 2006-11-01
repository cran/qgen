stat2e <- structure(function(paraDATA, alpha=0.05, v=list( c(-1,1),c(0,1),c(1,0),c(1,1),c(1,-1) ), Satterthwaite=TRUE, file=FALSE, path="~/Desktop/myproject.qgen"){
### default attributes          alpha=0.05; v=list( c(-1,1),c(0,1),c(1,0),c(1,1),c(1,-1) ); Satterthwaite=TRUE; file=FALSE
  ## CONDITION
  ## makes only sens for exactly two environments
  if (paraDATA@supl@enN==2 && !is.null(paraDATA)){
    ## reading
    part = paraDATA@orig@part
    hist = paraDATA@orig@hist
    chN = paraDATA@supl@chN
    enN = paraDATA@supl@enN
    ## fbN = paraDATA@supl@fbN
    ## rbN = paraDATA@supl@rbN
    siN = paraDATA@supl@siN
    daN = paraDATA@supl@daN
    idN = paraDATA@supl@idN
    fixe <- paraDATA@para@fixe # needed for the character and environment names
    ch.names <- dimnames(paraDATA@para@fixe)[[2]]
    en.names <- dimnames(paraDATA@para@fixe)[[1]]
    ##
    siS <- paraDATA@para@siS # sire variance-covariance matrix (upper case S for upper-case greek SIGMA)
    daS <- paraDATA@para@daS # dam variance-covariance matrix
    reS <- paraDATA@para@idS  +  diag(chN*enN) * paraDATA@para@error  +  2*siS - kronecker(diag(enN),matrix(1,nrow=chN, ncol=chN))*2*siS
    ## correspondaS to estimated phenotypic variance-covariance  matrix/ reSIGMA in the script
    rbS <- paraDATA@para@rbS # random block 
    phS <- paraDATA@para@phS # phenotypic (is only in the paraDATA if it is given)
    ## calculating the
    ##
    ## tblock
    ##    if( rbN>1){
    ## 	tSULd <- diag(tS)[1:chN]
    ## 	tSLRd <- diag(tS)[(chN+1):(2*chN)]
    ## 	tSLLd <- diag(matrix(tS[(chN+1):(2*chN) , (1:chN)],ncol=chN,nrow=chN)) # vector length=chN
    ## 	tSUL <- tS[1:chN,1:chN]
    ## 	tSLR <- tS[(chN+1):(2*chN),(chN+1):(2*chN)]
    ## 	tSLLLT <- tS[(chN+1):(2*chN),(1):(chN)]
    ## 	tSLLUT <- t(tS[(chN+1):(2*chN),(1):(chN)])
    ## 	tSVI <- tSULd+tSLRd-2*tSLLd #Variance of Inducibility
    ##     }
    ## sire
    siSULd <- diag(siS)[1:chN] #UpperLeftDiagonal ... Sigma environment 1
    siSLRd <- diag(siS)[(chN+1):(2*chN)] # Lower Right diagonal ... Sigma environment 2
    siSLLd <- diag(matrix(siS[(chN+1):(2*chN) , (1:chN)],ncol=chN,nrow=chN)) # vector length=chN...crosiS-environmental Covariance
    siSUL <- siS[1:chN, 1:chN] # ...Upper Left
    siSUR <- siS[1:chN, (chN+1):(2*chN)] # ... Upper Right
    siSLL <- siS[(chN+1):(2*chN), 1:chN] # ... Lower Left
    siSLR <- siS[(chN+1):(2*chN), (chN+1):(2*chN)] # ... Lower Right
    siSVI <- siSULd+siSLRd-2*siSLLd #Variance of Inducibility
    ## dam
    daSULd <- diag(daS)[1:chN]
    daSLRd <- diag(daS)[(chN+1):(2*chN)]
    daSLLd <- diag(matrix(daS[(chN+1):(2*chN) , (1:chN)],ncol=chN,nrow=chN)) # vector length=chN
    daSUL <- daS[1:chN,1:chN]
    daSUR <- daS[1:chN, (chN+1):(2*chN)] # ... Upper Right
    daSLL <- daS[(chN+1):(2*chN), 1:chN] # ... Lower Left
    daSLR <- daS[(chN+1):(2*chN),(chN+1):(2*chN)]
    daSVI <- daSULd+daSLRd-2*daSLLd #Variance of Inducibility
    ## individuals
    reSULd <- diag(reS)[1:chN]
    reSLRd <- diag(reS)[(chN+1):(2*chN)]
    reSLLd <- diag(matrix(reS[(chN+1):(2*chN) , (1:chN)],ncol=chN,nrow=chN)) # vector length=chN
    reSUL <- reS[1:chN,1:chN]
    reSUR <- reS[1:chN, (chN+1):(2*chN)] # ... Upper Right
    reSLL <- reS[(chN+1):(2*chN), 1:chN] # ... Lower Left
    reSLR <- reS[(chN+1):(2*chN),(chN+1):(2*chN)]
    reSVI <- reSULd+reSLRd-2*reSLLd #Variance of Inducibility
    ## Heritability of Plasticity
    herI <- 4*(siSVI) / (siSVI + daSVI + reSVI) # reS overritten if phenotype exists
### phenotypic; only theoretically!! NOT observable if every individual is measured in every environment
###
      if (!is.null(phS)) {
        phSULd <- diag(phS)[1:chN]
        phSLRd <- diag(phS)[(chN+1):(2*chN)]
        phSLLd <- diag(matrix(phS[(chN+1):(2*chN) , (1:chN)],ncol=chN,nrow=chN)) # vector length=chN
        phSUL <- phS[1:chN,1:chN]
        phSUR <- phS[1:chN, (chN+1):(2*chN)] # ... Upper Right
        phSLL <- phS[(chN+1):(2*chN), 1:chN] # ... Lower Left
        phSLR <- phS[(chN+1):(2*chN),(chN+1):(2*chN)]
        phSVI <- phSULd+phSLRd-2*phSLLd #Variance of Inducibility
        ## Heritability of Plasticity
        herI <- 4*(siSULd + siSLRd - 2*siSLLd) / phSVI
      }
### HERITABILITY with Confidence Intervals
### estimated phenotype
      herI <- matrix(herI, dimnames=list("Heritability of Inducibility",paste(ch.names,"",sep="")), nrow=1, ncol=chN)
      herI[herI>1] <- 0.999
### PERCENTAGE change
      if(hist[1]=="EMP") {
        fixx.percent <- NULL
        fix.utnf <- fixe
        for(cha in dimnames(fixe)[[2]]) {
          fix.utnf[,cha,] <- attr(dat.tra[,cha], "utnf")(fixe[,cha,])
          fixx.percent <- c(fixx.percent, mean((fix.utnf[1,cha,] - fix.utnf[2,cha,]) / fix.utnf[1,cha,]))
        }
        names(fixx.percent) <- dimnames(fixe)[[2]]
      }else{fix.utnf <- NULL;fixx.percent <- NULL}
### correlation inducibilities                                  -> corrI
      ScorrI <- (siSLR + siSUL - siSLL - siSUR) / sqrt(siSVI %*% t(siSVI))
      ScorrI <- matrix(ScorrI, dimnames=list(paste(ch.names,"INDUCIBILITY", sep="_"), paste(ch.names,"INDUCIBILITY", sep="_")), nrow=dim(ScorrI)[2])
      SpcorrI <- pcorr(ScorrI)
      PcorrI <- ((siSLR+daSLR+reSLR) + (siSUL+daSUL+reSUL) - (siSLL+daSLL+reSLL) - (siSUR+daSUR+reSUR)) / sqrt((siSVI+daSVI+reSVI) %*% t(siSVI+daSVI+reSVI))
      PcorrI <- matrix(PcorrI, dimnames=list(paste(ch.names,"INDUCIBILITY", sep="_"), paste(ch.names,"INDUCIBILITY", sep="_")), nrow=dim(PcorrI)[2])
      PpcorrI <- pcorr(PcorrI)
### correlation inducibility and one other trait                -> corrICa
      ScorrIC <- (siSLL-siSUL) / sqrt(siSVI %*% t(siSULd))
      ScorrIC <- matrix(ScorrIC, dimnames=list(paste(ch.names,"INDUCIBILITY", sep="_"), paste(en.names[1],ch.names,sep="_")), nrow=dim(ScorrIC)[2])
      ScorrIH <- (siSLR-siSUR) / sqrt(siSVI %*% t(siSLRd))
      ScorrIH <- matrix(ScorrIH, dimnames=list(paste(ch.names,"INDUCIBILITY", sep="_"), paste(en.names[2], ch.names,sep="_")), nrow=dim(ScorrIH)[2])
      PcorrIC <- ((siSLL+daSLL+reSLL)-(siSUL+daSUL+reSUL)) /sqrt((siSVI+daSVI+reSVI) %*% t(siSULd+daSULd+reSULd))
      PcorrIC <- matrix(PcorrIC, dimnames=list(paste(ch.names,"INDUCIBILITY", sep="_"), paste(ch.names,en.names[1],sep="_")), nrow=dim(ScorrIC)[2])
      PcorrIH <- ((siSLR+daSLR+reSLR)-(siSUR+daSUR+reSUR)) /sqrt((siSVI+daSVI+reSVI) %*% t(siSLRd+daSLRd+reSLRd))
      PcorrIH <- matrix(PcorrIH, dimnames=list(paste(ch.names,"INDUCIBILITY", sep="_"), paste(ch.names,en.names[2],sep="_")), nrow=dim(ScorrIH)[2])
      ScorrI[ScorrI>1] <- 0.9999
      SpcorrI[SpcorrI>1] <- 0.9999
      PcorrI[PcorrI>1] <- 0.9999
      PpcorrI[PpcorrI>1] <- 0.9999
      ScorrIC[ScorrIC>1] <- 0.9999
      ScorrIH[ScorrIH>1] <- 0.9999
      PcorrIC[PcorrIC>1] <- 0.9999
      PcorrIH[PcorrIH>1] <- 0.9999
      ScorrI[ScorrI<-1] <- -0.9999
      SpcorrI[SpcorrI<-1] <- -0.9999
      PcorrI[PcorrI<-1] <- -0.9999
      PpcorrI[PpcorrI < -1] <- -0.9999
      ScorrIC[ScorrIC < -1] <- -0.9999
      ScorrIH[ScorrIH < -1] <- -0.9999
      PcorrIC[PcorrIC < -1] <- -0.9999
      PcorrIH[PcorrIH < -1] <- -0.9999
###UNBAL
  if (length(paraDATA@spec@unbalanced)>0){ #+++unbal 
    ## CI after Bourdick92  p.107
    alpha <- 0.5*alpha
    w1u <- matrix(paraDATA@spec@unbalanced$w1u.vec,nrow=enN,ncol=chN) #dam in sireMS
    w2u <- matrix(paraDATA@spec@unbalanced$w2u.vec,nrow=enN,ncol=chN) #sire
    w3u <- matrix(paraDATA@spec@unbalanced$w3u.vec,nrow=enN,ncol=chN) #dam in damMS
    if(Satterthwaite){
      sDF <- matrix(paraDATA@spec@unbalanced$siDFappREML.vec,nrow=enN,ncol=chN)
      dDF <- matrix(paraDATA@spec@unbalanced$daDFappREML.vec,nrow=enN,ncol=chN)
    }else{
      sDF <- matrix(paraDATA@spec@unbalanced$siDF.vec,nrow=enN,ncol=chN)
      dDF <- matrix(paraDATA@spec@unbalanced$daDF.vec,nrow=enN,ncol=chN)
    } 
    iDF <- matrix(paraDATA@spec@unbalanced$idDF.vec,nrow=enN,ncol=chN)
#######
    ## The mean squares calculated with reml estimates
    ##aver <- function(x){mean(x)}
    aver <- function(x){2*x[1,]*x[2,]/(x[1,]+x[2,])}
    w1 <- aver(w1u); w2 <- aver(w2u); w3 <- aver(w3u) # aritmetic mean (to conservative ci)###555555
    S1U <- w2u[1,]*siSULd + w2u[2,]*siSLRd + w1u[1,]*daSULd + w1u[2,]*daSLRd + reSULd + reSLRd - 2*sqrt(w2u[1])*sqrt(w2u[2])*siSLLd - 2*sqrt(w1u[1])*sqrt(w1u[2])*daSLLd - 2*reSLLd  #; print(paste("S1U:", S1U))
    S2U <- w1u[1,]*daSULd + w1u[2,]*daSLRd + reSULd + reSLRd - 2*sqrt(w3u[1])*sqrt(w3u[2])*daSLLd - 2*reSLLd #; print(paste("S2U:", S2U))
    S3 <-  reSULd + reSLRd - 2*reSLLd   #; print(paste("S3:", S3))
    n1 <- (w2u[1,]*siSULd + w2u[2,]*siSLRd + w1u[1,]*daSULd + w1u[2,]*daSLRd + reSULd + reSLRd - 2*w2*siSLLd - 2*w1*daSLLd - 2*reSLLd)^2 / sum( (w2u[1,]^2*siSULd^2)/sDF[1,] + (w2u[2,]^2*siSLRd^2)/sDF[2,] + (w1u[1,]^2*daSULd^2)/dDF[1,] + (w1u[2,]^2*daSLRd^2)/dDF[2,] + reSULd^2/iDF[1,] + reSLRd^2/iDF[2,] )
    n2 <- (w3u[1,]*daSULd + w3u[2,]*daSLRd + reSULd + reSLRd - 2*w3*daSLLd - 2*reSLLd)^2 / sum( (w3u[1,]^2*daSULd^2)/dDF[1,] + (w3u[2,]^2*daSLRd^2)/dDF[2,] + reSULd^2/iDF[1,] + reSLRd^2/iDF[2,] )
    n3 <- (reSULd + reSLRd - 2*reSLLd)^2 / sum( reSULd^2/iDF[1,] + reSLRd^2/iDF[2,])
    ## CI after Bourdick92  p.107 (also the notation!); not from the original publication Sen92
    ## upper
      nomin <- w3*S1U - w1*qf(alpha,n1,n2)*S2U  -  (w3-w1)*qf(alpha,n1,n3)*S3
      denom <- w3*S1U - (w1-w2)*qf(alpha,n1,n2)*S2U  -  (w3-w1+w2-w2*w3)*qf(alpha,n1,n3)*S3
      CI.upper <- 4*nomin/denom
      ## lower
      nomin <- w3*S1U - w1*qf(1-alpha,n1,n2)*S2U  -  (w3-w1)*qf(1-alpha,n1,n3)*S3
      denom <- w3*S1U - (w1-w2)*qf(1-alpha,n1,n2)*S2U  -  (w3-w1+w2-w2*w3)*qf(1-alpha,n1,n3)*S3
      CI.lower <- 4*nomin/denom
      ## trimming
      for (oo in 1:chN){
        if (((w3*S1U)/(w1*S2U))[oo] < qf(1-alpha,n1[oo],n2[oo])) {CI.lower[oo] <- 0}
	if (((w3*S1U)/(w1*S2U))[oo] < qf(alpha,n1[oo],n2[oo])) {CI.upper[oo] <- 0}
        if(CI.lower[oo]>1){CI.lower[oo] <- 0.999999}
        if(CI.upper[oo]>1){CI.upper[oo] <- 0.999999}
        if(CI.upper[oo]<0){CI.upper[oo] <- 0}
      }
      CI.lower <- matrix(CI.lower, nrow=1, dimnames=list("Lower CI of plasticity", ch.names))
      CI.upper <- matrix(CI.upper, nrow=1, dimnames=list("Upper CI of plasticity", ch.names))
    }else{
      CI.upper <- NA; CI.lower <- NA
    }#___unbal
#      spec.out <- list(HER_ind=herI, fixx.percent=fixx.percent, fix.utnf=fix.utnf, Lower=CI.lower, Upper=CI.upper, sVI=siSVI, sC.i=ScorrI, sPC.i=SpcorrI, pC.i=PcorrI, pPC.i=PpcorrI, sC.ic=ScorrIC, sC.ih=ScorrIH, pC.ic=PcorrIC, pC.ih=PcorrIH) # old version
    }
### OUTPUT
###
  stat2e <- new("stat", orig=new("orig", hist=c(paraDATA@orig@hist, "stat2e"), warn=c(paraDATA@orig@warn, ""), time=c(paraDATA@orig@time, date())), stat=herI, lower.ci=CI.lower, upper.ci=CI.upper, lower.limes=0, upper.limes=1)
  ##
  if (file){
    save(stat2e, file=paste(path, "stat2e.rda", sep=""))
  }
  ##
  stat2
},stat.dim=1)
