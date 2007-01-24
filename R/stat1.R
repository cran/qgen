stat1 <- structure(function(paraDATA, alpha=0.05, frommethod="REML", Satterthwaite=TRUE, tex.table=FALSE, file=FALSE, path="~/qgen/") {
  ## frommethod (REML(default), ANOVA, ANOVAuw):
  ##      - only for experimental tests!
  ##      - whether the ANOVA estimators or
  ##        ANOVA estimators based on unweighted sums of squares should be used instead of REML
  ##      - only if awailable from the Est() function
  ##      - if the approximated degrees of freedom are NOT available,
  ##        the REML solutions are taken WITHOUT warning
  ## Satterthwaite (logical):
  ##      - only for experimental tests!
  ##      - whether the approximated sire and dam degree of freedom should be used
  ##      - if the approximated degrees of freedom are NOT available,
  ##        the standard degrees of freedom are taken WITHOUT warning
  ##
  ## attributes:           alpha=0.05; frommethod="REML";  Satterthwaite=TRUE; file=FALSE
### ALWAYS
###
  part = paraDATA@orig@part
  hist = paraDATA@orig@hist
  chN = paraDATA@supl@chN
  enN = paraDATA@supl@enN
  fbN = paraDATA@supl@fbN
  rbN = paraDATA@supl@rbN
  siN = paraDATA@supl@siN
  daN = paraDATA@supl@daN
  idN = paraDATA@supl@idN
  fixe <- paraDATA@para@fixe # only needed for character and environment names
###### POPULATION calculating the statistic from "known" population parameters
######
  if(identical(hist,"the") | identical(hist,c("the","sim")) ){#**** is POPulation or REAlization
    rbV <- diag(paraDATA@para@rbS)
    siV <- diag(paraDATA@para@siS)
    daV <- diag(paraDATA@para@daS)
    reV <- diag(paraDATA@para@idS)
    HER <- 4*siV/(rbV+siV+daV+reV)
#    HER <- matrix(HER, dimnames=list("Heritability", paste(rep(dimnames(fixe)[[1]],each=chN), "_",rep(dimnames(fixe)[[2]],times=enN), sep="")), nrow=1, ncol=length(HER))
    names(HER) <- paste(rep(dimnames(fixe)[[1]],each=chN), "_",rep(dimnames(fixe)[[2]],times=enN), sep="")
    ## dimnames(HER) <- list("Heritability",paste(rep(dimnames(fixx)[[1]],each=chN), "_",rep(dimnames(fixx)[[2]],times=enN), sep=""))
    ## trimming heritability
    HER[HER>0.999] <- 0.999
    HER[HER<0] <- 0
    ##output
    spec.out <- list(HER=HER, Lower=NULL, Upper=NULL)
    stat1 <- new("stat", orig=new("orig", hist=c(paraDATA@orig@hist, "stat1")), stat=HER, lower.ci=numeric(), upper.ci=numeric(), lower.limes=0, upper.limes=1)
  }else{# **** is POPulation or REAlization
###### ESTIMATION calculating the statistic from estimated parameters
######
### REML
###
  if(frommethod=="REML" &&(part=="REML" | part=="all") ){
    ## reading parData object
    rbV <- diag(paraDATA@para@rbS)
    siV <- diag(paraDATA@para@siS)
    daV <- diag(paraDATA@para@daS)
    ## for the reasoning behind the calculation of the residual variance see the documentation
    reV <- diag(paraDATA@para@idS  +  diag(chN*enN) * paraDATA@para@error  +  2*paraDATA@para@siS - kronecker(diag(enN),matrix(1,nrow=chN, ncol=chN))*2*paraDATA@para@siS) # corresponds to estimated phenotypic SIGMA/ rSIGMA in the script
    ## reading paraData for confidence Intervals
    ## testing for presence of ALPHA and UNBAL
    if (!is.null(alpha) && !is.null(paraDATA@spec@unbalanced) && rbV==0) {#### test for CI     IS THIS THE RIGHT TEST FOR THE EXISTENCE OF 
      ## Satterthwaite approximation
      if(Satterthwaite){
        n1 <- paraDATA@spec@unbalanced$siDFappREML.vec  # ; print(paste("n1appr",n1))
        n2 <- paraDATA@spec@unbalanced$daDFappREML.vec  # ; print(paste("n2appr",n2))
      }else{
        n1 <- paraDATA@spec@unbalanced$siDF.vec  # ; print(paste("n1",n1))
        n2 <- paraDATA@spec@unbalanced$daDF.vec  # ; print(paste("n2",n2))
      }
      n3 <- paraDATA@spec@unbalanced$idDF.vec
    } #___test for CI
  } #___REML
### ANOVA
###
  if (frommethod=="ANOVA"  && (part=="ANOVA" | part=="all") ){
    rbV <- zeroVector <- rep(0, times=chN*enN)
    siV <- paraDATA@spec$ANOVA@siS
    daV <- paraDATA@spec$ANOVA@daS
    reV <- paraDATA@spec$ANOVA@idS
    ##
     ## Satterthwaite approximation
    if(Satterthwaite){
      n1 <- paraDATA@spec@unbalanced$siDFappANOVA.vec  # ; print(paste("n1appr",n1))
      n2 <- paraDATA@spec@unbalanced$daDFappANOVA.vec  # ; print(paste("n2appr",n2))
    }else{
      n1 <- paraDATA@spec@unbalanced$siDF.vec  # ; print(paste("n1",n1))
      n2 <- paraDATA@spec@unbalanced$daDF.vec  # ; print(paste("n2",n2))
    }
    n3 <- paraDATA@spec@unbalanced$idDF.vec
  }
### ANOVAuw
###
  if (frommethod=="ANOVAuw"  && (part=="ANOVAuw" | part=="all") ){
    rbV <- zeroVector <- rep(0, times=chN*enN)
    siV <- paraDATA@spec$secondpartitioning$ANOVAuw$rand$siS
    daV <- paraDATA@spec$secondpartitioning$ANOVAuw$rand$daS
    reV <- paraDATA@spec$secondpartitioning$ANOVAuw$rand$idS
    ##
    ## Satterthwaite approximation
    if(Satterthwaite){
      n1 <- paraDATA@spec@unbalanced$siDFappANOVAuw.vec  # ; print(paste("n1appr",n1))
      n2 <- paraDATA@spec@unbalanced$daDFappANOVAuw.vec  # ; print(paste("n2appr",n2))
    }else{
      n1 <- paraDATA@spec@unbalanced$siDF.vec  # ; print(paste("n1",n1))
      n2 <- paraDATA@spec@unbalanced$daDF.vec  # ; print(paste("n2",n2))
    }
    n3 <- paraDATA@spec@unbalanced$idDF.vec
  }
### HERITABILITY calculation
###
  HER <- 4*siV/(rbV+siV+daV+reV)
  #HER <- matrix(HER, dimnames=list("Heritability",paste(rep(dimnames(fixe)[[1]],each=chN), "_",rep(dimnames(fixe)[[2]],times=enN), sep="")), nrow=1, ncol=length(HER))
  names(HER) <- paste(rep(dimnames(fixe)[[1]],each=chN), "_",rep(dimnames(fixe)[[2]],times=enN), sep="")
  ## dimnames(HER) <- list("Heritability",paste(rep(dimnames(fixx)[[1]],each=chN), "_",rep(dimnames(fixx)[[2]],times=enN), sep=""))
  ## trimming heritability
  HER[HER>0.999] <- 0.999
  HER[HER<0] <- 0
### CI calculation
###
  if (!is.null(alpha) && length(paraDATA@spec@unbalanced)!=0 && rbV==0) { #### test for CI
    ## reading 
    w1u <- paraDATA@spec@unbalanced$w1u.vec
    w2u <- paraDATA@spec@unbalanced$w2u.vec
    w3u <- paraDATA@spec@unbalanced$w3u.vec
    ##
    alpha <- 0.5*alpha
    ## The mean squares calculated with reml estimates
    S1U <- reV + w1u*daV + w2u*siV   #; print(paste("S1U:", S1U))
    S2U <- reV + w3u*daV   #; print(paste("S2U:", S2U))
    S3 <- reV   #; print(paste("S3:", S3))
    ## CI after Bourdick92  p.107 (also the notation!); not from the original publication Sen92
    ## upper
    nomin <- w3u*S1U - w1u*qf(alpha,n1,n2)*S2U  -  (w3u-w1u)*qf(alpha,n1,n3)*S3   #;  print(paste("1:", w3u*S1U, "2:", w1u*qf(alpha,n1,n2)*S2U,  "3:",(w3u-w1u)*qf(alpha,n1,n3)*S3))
    denom <- w3u*S1U - (w1u-w2u)*qf(alpha,n1,n2)*S2U  -  (w3u-w1u+w2u-w2u*w3u)*qf(alpha,n1,n3)*S3#;  print(paste("1:", w3u*S1U, "2:", (w1u-w2u)*qf(alpha,n1,n2)*S2U,  "3:",(w3u-w1u+w2u-w2u*w3u)*qf(alpha,n1,n3)*S3))
    remlU <- 4*nomin/denom
    ## lower
    nomin <- w3u*S1U - w1u*qf(1-alpha,n1,n2)*S2U  -  (w3u-w1u)*qf(1-alpha,n1,n3)*S3#;  print(paste("1:", w3u*S1U, "2:", w1u*qf(alpha,n1,n2)*S2U,  "3:",(w3u-w1u)*qf(alpha,n1,n3)*S3))
    denom <- w3u*S1U - (w1u-w2u)*qf(1-alpha,n1,n2)*S2U  -  (w3u-w1u+w2u-w2u*w3u)*qf(1-alpha,n1,n3)*S3#;  print(paste("1:", w3u*S1U, "2:", (w1u-w2u)*qf(alpha,n1,n2)*S2U,  "3:",(w3u-w1u+w2u-w2u*w3u)*qf(alpha,n1,n3)*S3))
    remlL <- 4*nomin/denom
    ## trimming
    for (oo in 1:length(n1)) {
      if ( ((w3u*S1U) / (w1u*S2U))[oo] < qf(1-alpha,n1[oo],n2[oo]) ) {remlL[oo] <- 0}
      if ( ((w3u*S1U) / (w1u*S2U))[oo] < qf(alpha,n1[oo],n2[oo]) ) {remlU[oo] <- 0}
      if(remlL[oo]>1){remlL[oo] <- 1}
      if(remlU[oo]>1){remlU[oo] <- 1}
      if(remlL[oo]<0){remlL[oo] <- 0}
      if(remlU[oo]<0){remlU[oo] <- 0}
    }
### Burdick p.101
    G1=1-1/qf(1-alpha,n1,Inf)
    G2=1-1/qf(1-alpha,n2,Inf)
    G3=1-1/qf(1-alpha,n3,Inf)
    H1=1/qf(alpha,n1,Inf) -1
    H2=1/qf(alpha,n2,Inf) -1
    H3=1/qf(alpha,n3,Inf) -1
    G12=((qf(1-alpha,n1,n2)-1)^2  -  G1^2*(qf(1-alpha,n1,n2))^2 - H2^2)   /    (qf(1-alpha,n1,n2))
    G13=((qf(1-alpha,n1,n3)-1)^2  -  G1^2*(qf(1-alpha,n1,n3))^2 - H3^2)   /    (qf(1-alpha,n1,n3))
    G32=((qf(1-alpha,n3,n2)-1)^2  -  G3^2*(qf(1-alpha,n3,n2))^2 - H2^2)   /    (qf(1-alpha,n3,n2))
    G13star=(1-1/qf(1-alpha,n1+n3,Inf))^2 * (n1+n3)^2/(n1*n3) - G1^2*n1/n3 - G3^2*n3/n1
    H12= ((1-qf(alpha, n1, n2))^2 - H1^2*qf(alpha,n1,n2)^2 - G2^2)  /    (qf(alpha, n1, n2))
    H13= ((1-qf(alpha, n1, n3))^2 - H1^2*qf(alpha,n1,n3)^2 - G3^2)  /    (qf(alpha, n1, n3))
    H32= ((1-qf(alpha, n3, n2))^2 - H3^2*qf(alpha,n3,n2)^2 - G2^2)  /    (qf(alpha, n3, n2))
    H23star=(1-1/qf(1-alpha,n2+n3,Inf))^2 * (n2+n3)^2/(n2*n3) - G2^2*n2/n3 - G3^2*n3/n1
    c2=w1u/w3u
    c3=c2-1
    VL <- G1^2*S1U^2 + G3^2*c3*S3^2 + H2^2*c2^2*S2U^2 + G12*c2*S1U*S2U + G32*c3*c2*S3*S2U + G13star*c3*S1U*S3
    VL[c3<0] <- G1[c3<0]^2*S1U[c3<0]^2 + H2[c3<0]^2*c2[c3<0]*S2U[c3<0]^2 + H3[c3<0]^2*c3[c3<0]^2*S3[c3<0]^2 + G12[c3<0]*c2[c3<0]*S1U[c3<0]*S2U[c3<0] + G13[c3<0]*abs(c3)[c3<0]*S1U[c3<0]*S3[c3<0]
    VU <- H1^2*S1U^2 + H3^2*c3*S3^2 + G2^2*c2^2*S2U^2 + H12*c2*S1U*S2U + H32*c3*c2*S3*S2U
    VU[c3<0] <- H1[c3<0]^2*S1U[c3<0]^2 + G2[c3<0]^2*c2[c3<0]*S2U[c3<0]^2 + G3[c3<0]^2*c3[c3<0]^2*S3[c3<0]^2 + H12[c3<0]*c2[c3<0]*S1U[c3<0]*S2U[c3<0] + H13[c3<0]*abs(c3)[c3<0]*S1U[c3<0]*S3[c3<0] + H23star[c3<0]*abs(c3)[c3<0]*S2U[c3<0]*S3[c3<0]
    ##
    CIsiV.lower <- (S1U - c2*S2U + c3*S3 - sqrt(VL))   /  (w2u)
    CIsiV.upper <- (S1U - c2*S2U + c3*S3 + sqrt(VU))   /  (w2u)
    ##trimming
    for (oo in 1:length(n1)) {
      if(CIsiV.lower[oo]>1){CIsiV.lower[oo] <- 1}
      if(CIsiV.upper[oo]>1){CIsiV.upper[oo] <- 1}
      if(CIsiV.lower[oo]<0){CIsiV.lower[oo] <- 0}
      if(CIsiV.upper[oo]<0){CIsiV.upper[oo] <- 0}
      if(tex.table){
        cat(paste("\\ci{",format(round(remlL[oo],3),nsmall=2),"}{",format(round(remlU[oo],3),nsmall=2),"}\n"), file=paste(path,"CIhertex.tex",sep=""), append=TRUE)
        cat(paste("&","\n","\\ci{",format(round(CIsiV.lower[oo],6),nsmall=3),"}{",format(round(CIsiV.upper[oo],3),nsmall=3),"}\n"), file=paste(path,"CIVtex.tex", sep=""), append=TRUE)
      }
    }
#     ## this part is just for the curious!
#     ## for balanced data this should give the same resuts as the code just above
#     ## after Bourdick92  p.86 (original: Graybill and Wang 1979)
#     remlU <- NULL
#     remlL <- NULL
#     alpha <- 0.975
#     for (oo in 1:length(n1)){
#       print(S1U[oo]/S2U[oo])
#       print(qf(alpha,n1[oo],n2[oo]))
#       if (S1U[oo]/S2U[oo] <= qf(alpha,n1[oo],n2[oo])){remlL.oo <- 0} else {
#         remlL.oo <- 4 * (S1U[oo] - qf(alpha,n1[oo],n2[oo])*S2U[oo])    /     (S1U[oo] + (design[5]-1)*qf(alpha,n1[oo],n2[oo])*S2U[oo] + design[5]*(design[6]-1)*qf(alpha,n1[oo],n3[oo])*S3[oo])
#       }
#       alpha <- 1-alpha
#       remlU.oo <- 4 * (S1U[oo] - qf(alpha,n1[oo],n2[oo])*S2U[oo])    /     (S1U[oo] + (design[5]-1)*qf(alpha,n1[oo],n2[oo])*S2U[oo] + design[5]*(design[6]-1)*qf(alpha,n1[oo],n3[oo])*S3[oo])
#       remlU <- c(remlU,remlU.oo)
#       remlL <- c(remlL,remlL.oo)
#     }
    ##output
  #  spec.out <- list(HER=HER, Lower=remlL, Upper=remlU, CIsiV.lower=CIsiV.lower, CIsiV.upper=CIsiV.upper)
    stat1 <- new("stat", orig=new("orig", hist=c(paraDATA@orig@hist, "stat1"), warn=c(paraDATA@orig@warn, ""), time=c(paraDATA@orig@time, date())), stat=HER, lower.ci=remlL, upper.ci=remlU, lower.limes=0, upper.limes=1)
  }else{
    stat1 <- new("stat", orig=new("orig", hist=c(paraDATA@orig@hist, "stat1"), warn=c(paraDATA@orig@warn, ""), time=c(paraDATA@orig@time, date())), stat=HER, lower.ci=numeric(), upper.ci=numeric(), lower.limes=0, upper.limes=1)
    #spec.out <- list(HER=HER, Lower=NULL, Upper=NULL)
  } #___ test for CI
### untransformation of the Coefficitents
### calculates by how many percentages the character value decreased from environment 1 to 2
### (on the original scale)
#   if(any(paraDATA$orig$hist=="Emp")&enN==2) {
#     fixx <- paraDATA$para$fix
#     fixx.percent <- NULL
#     fix.utnf <- fixx
#     for(cha in dimnames(fixx)[[2]]) {
#       fix.utnf[,cha,] <- attr(dat.tra[,cha], "utnf")(fixx[,cha,])
#       fixx.percent <- c(fixx.percent, mean((fix.utnf[1,cha,] - fix.utnf[2,cha,]) / fix.utnf[1,cha,]))
#     }
#     ##output
#     perc.out <- list(fix.utnf=fix.utnf, fix.percent=fixx.percent)
#   }else{fix.utnf <- NULL;fixx.percent <- NULL}
} #___POPulation

### OUTPUT
###
  ## orig
  #origStat1 <- new("orig", hist=c(paraDATA@orig@hist, "Stat1"))
  ## stat
  #stat1 <- new("stat", orig=origStat1, stat=HER, lower.ci=numeric(), upper.ci=numeric(), lower.limes=0, upper.limes=1)
#   ## supl
#   suplStat1 <- paraDATA@supl
#   para <- NULL
#   data <- NULL
#   spec <- spec.out
#   ##
#   stat1 <- list(orig=orig, supl=supl, para=para, data=data, spec=spec)
  ##
  if (file){
    dir.create(path=path, showWarnings = FALSE, recursive = TRUE)
    save(stat1, file=paste(path, "stat1.rda", sep=""))
  }
  stat1
},stat.dim=1)
########################################################################
### Sen92
#   alpha <- 0.5*alpha # because alpha was two-sided
#   S <- paraDATA$meth@unbalanced$S.vec
#   D.s <- paraDATA$meth@unbalanced$D.s.vec
# #  I.s <- paraDATA$meth@unbalanced$I.s.vec
#   I.sd <- paraDATA$meth@unbalanced$I.sd.vec
#   H.s <- paraDATA$meth@unbalanced$HP.vec
### Numerical Example Sen92;p272
#   S <- 4
#   D.s <- c(2,3,4,5)
#   I.sd <- c(2,3,2,3,4,2,3,4,5,2,2,3,4,5)
  #
#   k1 <- S^-1 * sum(D.s^-1)
#   k2 <- S^-1 * sum( (D.s^-1)*(H.s^-1) )
#   k3 <- sum( (D.s-1)/(sum(D.s)*H.s) )
#   SA <- rV*k2 + dV*k1 + sV
#   SB <- rV*k3 + dV
#   SC <- rV
#   sigmaA <- max(0,SA - SB*1/S*sum(D.s^-1) - SC* (1/S*sum( (D.s^-1)*(H.s^-1) ) - sum( (D.s-1)/(sum(D.s)-S) *H.s^-1 ) * ( (1/S)*sum(D.s^-1))))
#   sigmaB <- max(0, SB - SC*sum((D.s-1)/(sum(D.s)*H.s)))
#   sigmaC <- SC
# #   SA <- rV + 2*dV + 4*sV
# #   SB <- rV + 2*dV
# #   SC <- rV
# #   sigmaA <- sV
# #   sigmaB <- dV
# #   sigmaC <- rV
#   ap <- D.s^-1 - (sum(D.s^-1)/S)
#   bp <- D.s^-1*H.s^-1 - sum(D.s^-1*H.s^-1)/S
#      H.sd <- rep(H.s, times=D.s) # it has to have the full length
#      D.sd <- rep(D.s, times=D.s) # it has to have the full length
#   c1 <- sum(( (D.sd-2)/(sum(D.s)-S) ) * ( ((I.sd^-1) * (H.sd^-1))^2 /D.sd ) )
#   c2 <- sum( (D.s-1)/(sum(D.s)-S) * ( (H.s^-1) - sum((D.s-1)/(sum(D.s)*H.s)) )^2)
#   nA <- (S-1)*(1 + S*(S-2)/S-1*sum((ap*sigmaB+bp*sigmaC)^2)/(sum(sigmaA + D.s^-1*sigmaB + D.s^-1*H.s^-1*sigmaC))^2)^-1
#   nB <- (sum(D.s)-S)*(1 + sigmaC^2*(c1+c2)/sigmaB^2)^-1
#   nC <- sum(I.sd) - sum(D.s)
#   nomin <- SA - SB*qf(alpha,nB,nA)*k1 - SC*qf(alpha,nC,nA)*(k2 - k3*k1)
#   denom <- SA- SB*qf(alpha,nB,nA)*(k1-1) - SC*(k2-k3*k1+k3-1)*qf(alpha,nC,nA)
#   remlL <- NULL
#   for (oo in 1:length(nA)) {
#     if ( (SA/SB) > qf(alpha,nB[oo],nA[oo]) ) {remlL[oo] <- 0} else {remlL[oo] <- 4*nomin/denom}
#   }
#   alpha <- 1-alpha
#   nomin <- SA - SB*qf(alpha,nB,nA)*k1 - SC*qf(alpha,nC,nA)*(k2 - k3*k1)
#   denom <- SA- SB*qf(alpha,nB,nA)*(k1-1) - SC*(k2-k3*k1+k3-1)*qf(alpha,nC,nA)
#   remlU <- NULL
#   for (oo in 1:length(nA)) {
#     if ( (SA/SB) > qf(alpha,nB[oo],nA[oo]) ) {remlU[oo] <- 0} else {remlU[oo] <- 4*nomin/denom}
#   }
