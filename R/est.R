est <- function(paraDATA, file=TRUE, path="~/Desktop/myproject.qgen/"){
##################################################
###                    1 or 2 environment VERSIONS
  ## under the assumption that with a newer Version of R also the packages were updated!
  if(R.Version()$minor>2){
    versionOLD <- TRUE
  }
  versionOLD <- FALSE
  print(paste("is this the old version?", as.character(versionOLD)))
##################################################
###                                        reading
  chN <- paraDATA@supl@chN
  enN <- paraDATA@supl@enN
  fbN <- paraDATA@supl@fbN
  rbN <- paraDATA@supl@rbN
  siN <- paraDATA@supl@siN
  daN <- paraDATA@supl@daN
  idN <- paraDATA@supl@idN
  fixe <- paraDATA@supl@miss #for names!
  dat <- paraDATA@DATA@dat
  dat$si_da <- as.factor(paste(dat$si,"_", dat$da, sep=""))
  dat$en_ch<- as.factor(paste(dat$en,"_", dat$ch, sep=""))
##################################################
  ## deciding what should be calculated, depending on the origin (part and hist)
  part <- paraDATA@orig@part
  if(part=="all"){
    REMLpartitioning <- TRUE
    ANOVApartitioning <- TRUE
  }
  if(part=="REML"){
    REMLpartitioning <- TRUE
    ANOVApartitioning <- FALSE
  }
  if(chN*enN==1 &&(part=="ANOVA" | part=="ANOVAuw")){
    ## ANOVA is also possible for chN*enN>1 but ONLY as second partitioning
    ## then for every combination seperately
    REMLpartitioning <- FALSE
    ANOVApartitioning <- TRUE
  }else{
    REMLpartitioning <- TRUE
    ANOVApartitioning <- FALSE
  }
  ##
  hist <- paraDATA@orig@hist
  if(identical(hist, "emp")){ #"full"
    ANOVApartitioning <- TRUE
    unbalanced <- TRUE
    secondcontrast <- ifelse(enN==2, TRUE, FALSE) # because secondcontrast only makes sense to compare twoenvironments
    modelsummary <- TRUE
  }
  if(identical(hist, c("the", "sim"))){ #"S"
    ANOVApartitioning <- TRUE
    unbalanced <- TRUE
    secondcontrast <- FALSE
    modelsummary <- FALSE
  }
  if(identical(hist, c("emp", "est", "sim")) | identical(hist, c("the", "sim", "est", "sim")) ){  #"R"
    unbalanced <- FALSE
    secondcontrast <- FALSE
    modelsummary <- FALSE
  }
  if(identical(hist, c("emp", "est", "sim", "est", "sim")) | identical(hist, c("the", "sim", "est", "sim", "est", "sim"))){ #"Q"
    unbalanced <- FALSE
    secondcontrast <- FALSE
    modelsummary <- FALSE
  }
##################################################
###                                           REML
### always if(unbalanced) also if(ANOVApartitioning)
  if(REMLpartitioning){#(REML+++)
    ##
    if(!versionOLD){#(+++newversion) # lme4version0.995-2
      if(chN*enN*fbN ==1){#++1a
        model <- lmer(y ~ 1 + (1|si) + (1|si_da), dat)
        sBLUP <- ranef(model)$si # SIRE BLUPs
        names(sBLUP) <- names(sBLUP)
        sBLUP <- data.matrix(sBLUP)
        fix.coef <- fixe
      }else{#--1a
        stop("This version of qgen does not support calculations with more than one environment, character or fixed block!")
        dat$fbchen = factor(paste(dat$fb,dat$ch,dat$en, sep="_")) # defining a new factor
        model <- lmer(y ~ fbchen + (en_ch-1|si) + (en_ch-1|si_da) + (en_ch-1|id), dat) # THE HEART !!
        fix.coef <- array(fixef(model)[1] + contrasts(dat$fbchen)%*% fixef(model)[-1], dim=c(length(levels(dat$en)), length(levels(dat$ch)), length(levels(dat$fb))), dimnames=list(levels(dat$en), levels(dat$ch), levels(dat$fb)))
      }#__1a
      ## getting all parameters from the random part of the lme4object ("model")
      scale <- attributes(VarCorr(model))$sc
      ## sire
      sireSigma <- as.matrix(VarCorr(model)$si)
      n <- length(sireSigma) # used later to form matrices
      sireSigma <- matrix(sireSigma, dimnames=list(paste(rep(levels(dat$en),each=chN),"_",rep(levels(dat$ch),times=enN),sep=""),paste(rep(levels(dat$en),each=chN),"_",rep(levels(dat$ch),times=enN),sep="") ), ncol=n, nrow=n)
      ## dam
      damSigma <- as.matrix(VarCorr(model)$si_da)
      damSigma <- matrix(damSigma, dimnames=list(paste(rep(levels(dat$en),each=chN),"_",rep(levels(dat$ch),times=enN),sep=""),paste(rep(levels(dat$en),each=chN),"_",rep(levels(dat$ch),times=enN),sep="") ), ncol=n, nrow=n)
      ## individual
      if ((chN*enN)>1){
        iSD <- VarCorr(model)$id@x
        indSigma <- matrix(iSD,ncol=n,nrow=n,byrow=FALSE) * matrix(iSD,ncol=n,nrow=n,byrow=TRUE) * VarCorr(model)@reSumry$id@.Data
        indSigma <- matrix(indSigma, dimnames=list(paste(rep(levels(dat$en),each=chN),"_",rep(levels(dat$ch),times=enN),sep=""),paste(rep(levels(dat$en),each=chN),"_",rep(levels(dat$ch),times=enN),sep="") ), ncol=n, nrow=n)
      }else{
        indSigma <- matrix(0,nrow=n,ncol=n)
        dimnames(indSigma) <- list(paste(rep(levels(dat$en),each=chN),"_",rep(levels(dat$ch),times=enN),sep=""),paste(rep(levels(dat$en),each=chN),"_",rep(levels(dat$ch),times=enN),sep="") )
        eS <- scale^2
        iSD <- 0
        eC <- NULL
        ePC <- NULL
      }
      ## different blocking factors are not yet supported
      if(3==4){ #(rbN>1){
        tSD <- VarCorr(model)@reSumry$rb@stdDev * scale
        tblockSigma <- matrix(tSD,ncol=n,nrow=n,byrow=FALSE) * matrix(tSD,ncol=n,nrow=n,byrow=TRUE) * VarCorr(model)@reSumry$rb@.Data
        tblockSigma <- matrix(tblockSigma, dimnames=list(paste(rep(levels(dat$en),each=chN),"_",rep(levels(dat$ch),times=enN),sep=""),paste(rep(levels(dat$en),each=chN),"_",rep(levels(dat$ch),times=enN),sep="") ), ncol=n, nrow=n)
      }else{
        tblockSigma <- matrix(0, dimnames=list(paste(rep(levels(dat$en),each=chN),"_",rep(levels(dat$ch),times=enN),sep=""),paste(rep(levels(dat$en),each=chN),"_",rep(levels(dat$ch),times=enN),sep="") ), ncol=n,nrow=n)
        tC <- NULL
        tPC <- NULL
      }
    }#(___newversion)
### ---------------------------------
    if(versionOLD){# (+++oldversion)
      ## lme4version0.95-6
      if(chN*enN*fbN ==1){#(1env+++)
        model <- lmer(y ~ 1 + (1|si) + (1|si_da), dat)
        summary(model)
        sBLUP <- ranef(model)$si # SIRE BLUPs
        names(sBLUP) <- names(sBLUP)
        sBLUP <- data.matrix(sBLUP)
        fix.coef <- fixe
      }else{#(1env---)
        dat$fbchen = factor(paste(dat$fb,dat$ch,dat$en, sep="_")) # defining a new factor
        model <- lmer(y ~ fbchen + (en_ch-1|si) + (en_ch-1|si_da) + (en_ch-1|id), dat) # THE HEART !!
        summary(model)
        fix.coef <- array(fixef(model)[1] + contrasts(dat$fbchen)%*% fixef(model)[-1], dim=c(length(levels(dat$en)), length(levels(dat$ch)), length(levels(dat$fb))), dimnames=list(levels(dat$en), levels(dat$ch), levels(dat$fb)))
      }#(1env___)
      scale <- VarCorr(model)@scale
      sSD <- VarCorr(model)@reSumry$si@stdDev * scale
      dSD <- VarCorr(model)@reSumry$si_da@stdDev * scale
      n <- length(sSD)
      sireSigma <- matrix(sSD,ncol=n,nrow=n,byrow=FALSE) * matrix(sSD,ncol=n,nrow=n,byrow=TRUE) * VarCorr(model)@reSumry$si@.Data
      sireSigma <- matrix(sireSigma, dimnames=list(paste(rep(levels(dat$en),each=chN),"_",rep(levels(dat$ch),times=enN),sep=""),paste(rep(levels(dat$en),each=chN),"_",rep(levels(dat$ch),times=enN),sep="") ), ncol=n, nrow=n)
      damSigma <- matrix(dSD,ncol=n,nrow=n,byrow=FALSE) * matrix(dSD,ncol=n,nrow=n,byrow=TRUE) * VarCorr(model)@reSumry$si_da@.Data
      damSigma <- matrix(damSigma, dimnames=list(paste(rep(levels(dat$en),each=chN),"_",rep(levels(dat$ch),times=enN),sep=""),paste(rep(levels(dat$en),each=chN),"_",rep(levels(dat$ch),times=enN),sep="") ), ncol=n, nrow=n)
      if ((chN*enN)>1){ # id
        iSD <- VarCorr(model)@reSumry$id@stdDev * scale
        indSigma <- matrix(iSD,ncol=n,nrow=n,byrow=FALSE) * matrix(iSD,ncol=n,nrow=n,byrow=TRUE) * VarCorr(model)@reSumry$id@.Data
        indSigma <- matrix(indSigma, dimnames=list(paste(rep(levels(dat$en),each=chN),"_",rep(levels(dat$ch),times=enN),sep=""),paste(rep(levels(dat$en),each=chN),"_",rep(levels(dat$ch),times=enN),sep="") ), ncol=n, nrow=n)
      }else{
        indSigma <- matrix(0,nrow=n,ncol=n)
        dimnames(indSigma) <- list(paste(rep(levels(dat$en),each=chN),"_",rep(levels(dat$ch),times=enN),sep=""),paste(rep(levels(dat$en),each=chN),"_",rep(levels(dat$ch),times=enN),sep="") )
        eS <- scale^2
        iSD <- 0
        eC <- NULL
        ePC <- NULL
      }
      if(3==4){ #(rbN>1){# tb
        tSD <- VarCorr(model)@reSumry$rb@stdDev * scale
        tblockSigma <- matrix(tSD,ncol=n,nrow=n,byrow=FALSE) * matrix(tSD,ncol=n,nrow=n,byrow=TRUE) * VarCorr(model)@reSumry$rb@.Data
        tblockSigma <- matrix(tblockSigma, dimnames=list(paste(rep(levels(dat$en),each=chN),"_",rep(levels(dat$ch),times=enN),sep=""),paste(rep(levels(dat$en),each=chN),"_",rep(levels(dat$ch),times=enN),sep="") ), ncol=n, nrow=n)
      }else{
        tblockSigma <- matrix(0, dimnames=list(paste(rep(levels(dat$en),each=chN),"_",rep(levels(dat$ch),times=enN),sep=""),paste(rep(levels(dat$en),each=chN),"_",rep(levels(dat$ch),times=enN),sep="") ), ncol=n,nrow=n); tC <- NULL; tPC <- NULL
      }
    }#(___oldversion)
    ##output
    reml.para <- new("para", rbS=tblockSigma, siS=sireSigma, daS=damSigma, idS=indSigma, phS=matrix(), error=scale^2, fixe=fix.coef)
#################################
### LONGOUTPUT
    if(modelsummary){
      sBLUP <- ranef(model)$si # SIRE BLUPs (best linear unbiased predictors)
      names(sBLUP) <- substr(names(sBLUP),5,99)
      sBLUP <- data.matrix(sBLUP) # a matrix with sires in rows and characters in columns
      ##output
    modelsummary.reml <- list(sBLUP=sBLUP, modelsumry=summary(model))
    }else{
      modelsummary.reml <- list()
    }
#################################
### SECONDCONTRAST
    ## a second model with a different contrast matrix
    if (secondcontrast){#(secondcontrast+++)
      dat$chenfb <- factor(paste(dat$ch,dat$en,dat$fb, sep="_")) # so we can use the contr.sdif from MASS to get a test in every trait of the environmental effect (as long as there are only 2 environments)
      # DIESE Aenderung nochmals testen!! neu mit dat$fbchen
      ## CONTRASTS
      ntb <- length(levels(dat$fb))
      nen <- length(levels(dat$en))
      ntr <- length(levels(dat$ch))
      ## contrasts for environmental differences
      i <- diag(1, nrow=ntr)
      e <- rep(c(-1,1),each=ntb)
      k.endiff <- kronecker(i,e)
      dimnames(k.endiff) <- list(NULL,paste(levels(dat$ch),"_enDIFF",sep=""))
      ## contrasts for trait deviance from overall mean
      mm <- matrix(-1, ntr, ntr)
      diag(mm) <-ntr-1
      k.mean <- as.matrix(kronecker(mm,rep(1,times=ntb*nen))[,-1])
      dimnames(k.mean) <- list(NULL,levels(dat$ch)[-1])
      ## contrasts for differences among tb's; within trait and environment
      if(ntb>1){
        itb <- diag(1, nrow=ntr*nen)
        etb <- matrix(rep(c(c(1,-1),rep(0,times=ntb-1)),ntb-1)[1:(ntb*(ntb-1))], nrow=ntb, ncol=ntb-1)
        k.tbdiff <- kronecker(itb,etb)
        dimnames(k.tbdiff) <- list(NULL,as.character(1:(dim(k.tbdiff)[2])))
      }else{
        k.tbdiff <- NULL
      }
      ##
      k <- cbind(k.mean, k.endiff, k.tbdiff)[,1:(dim(k.endiff)[1]-1)]
      contrasts(dat$chenfb) <- k
      model2 <- lmer(y ~ chenfb + (en_ch-1|si) + (en-1|si_da) + (en_ch-1|id), dat) # THE second HEART !!
      ## a matirx with all environment - trait combinations (mean over all tblocks)
      fix.coef2 <- fixef(model2) 
      int <- fix.coef2[1]
      if (ntr >1){tr.means <- c(int - fix.coef2[2], fix.coef2[2:(1+ntr-1)]+int)}else{tr.means <- c(int-fix.coef2[2])}
      names(tr.means) <- levels(dat$ch)
      en.diff <- fix.coef2[(1+ntr):(2*ntr)]
      names(en.diff) <- levels(dat$ch)
      FIX <- rbind(tr.means-en.diff, tr.means+en.diff) # a matrix with environments as rows and characters as columns
      dimnames(FIX) <- list(levels(dat$en),levels(dat$ch)) #here the names are readded ATTENTION!
      ##output
      secondcontrast.reml <- list(FIX=FIX, modelsumry2=model2)
    }else{
      secondcontrast.reml <- list() 
    }#(secondcontrast___)
  }#(REML___)
### UNBALANCED
### 
  if(unbalanced | ANOVApartitioning) {#(unbalanced+++)
    dat$sidaid <- factor(paste(dat$si,"_",dat$da,"_",dat$id,sep=""),ordered=TRUE)
    siDF.vec <- NULL; daDF.vec <- NULL; idDF.vec <- NULL
    siDFappREML.vec <- NULL; daDFappREML.vec <- NULL
    siDFappANOVA.vec <- NULL; daDFappANOVA.vec <- NULL
    siDFappANOVAuw.vec <- NULL; daDFappANOVAuw.vec <- NULL
    w1u.vec <- NULL; w2u.vec <- NULL; w3u.vec <- NULL
    reV.ANOVA.vec <- NULL; daV.ANOVA.vec <- NULL; siV.ANOVA.vec <- NULL
    reV.ANOVAuw.vec <- NULL; daV.ANOVAuw.vec <- NULL; siV.ANOVAuw.vec <- NULL
    ##for all trt-trait combinations
    for (i in levels(dat$en_ch)) { # e.g.   i <- "en1_ch1"
      ##  
      dat.tren <- dat[dat$en_ch==i,]
      dat.tren$si <- factor(dat.tren$si, exclude=NULL) # if not all sires are represented we need to adjust the number of levels
      dat.tren$si <- factor(paste(dat.tren$si,"_","si",sep="")) # this sorts them "correctly" like sidain
      dat.tren$da <- factor(paste(dat.tren$da,"_","da",sep="")) # this sorts them "correctly" like sidain
      dat.tren$id <- factor(paste(dat.tren$id,"_","id",sep="")) # this sorts them "correctly" like sidain
      ##Numerical Example Sen92;p272
      ## S
      S <- length(levels(dat.tren$si)) #total number of sires
      ## I.sd
      n.ij <- tapply(dat.tren$id, dat.tren$si_da, function(x) length(tabulate(x)[tabulate(x)>0]))[order(levels(dat.tren$sida))]  # number of ind within sire and dam
      I.sd <- n.ij[!is.na(n.ij)]
      ## I.s
      n.i <- tapply(dat.tren$id, dat.tren$si, function(x) length(tabulate(x)[tabulate(x)>0]))[order(levels(dat.tren$si))]
      I.s <- n.i[!is.na(n.i)]   # number of ind within sire
      ## D.s
      M.i <- tapply(dat.tren$da, dat.tren$si, function(x) length(tabulate(x)[tabulate(x)>0])) #number of dams within sire
      M.i <- M.i[!is.na(M.i)]
      D.s <- M.i[order(levels(dat.tren$si))] # number of dams within sire
      D.s <- D.s[!is.na(D.s)]
      ##
      H.s <- D.s/tapply(1/I.sd, rep(sort(levels(dat.tren$si)), times=D.s),sum)  #2.1
      ## note:      rep(sort(levels(dat.tren$si)), times=D.s) # repeats the sires to the total length of sida
      #
      siDF <- S-1
      daDF <- sum(D.s)-S##
      idDF <- sum(I.sd) - sum(D.s)
      w1u <- sum(1/D.s) / sum(1/(D.s*H.s))
      w2u <- S/sum(1/(D.s*H.s))
      w3u <- daDF/sum((D.s-1)/H.s)
       ## to be exported
      siDF.vec <- c(siDF.vec, siDF)
      daDF.vec <- c(daDF.vec, daDF)
      idDF.vec <- c(idDF.vec, idDF)
      w1u.vec <- c(w1u.vec, w1u)
      w2u.vec <- c(w2u.vec, w2u)
      w3u.vec <- c(w3u.vec, w3u)
      ## Satterthwaite approximation (see Sen1992; equations 2.5 and 2.6; NOT described in Burdick 1992)
      ## to estimate the degree of freedom for the factors sire and dam (nA, nB), which are exported and then used (e.g. in Stat1())
      ## the estimated variance components are needed to calculate the Satterthwaite approximation of the sire and dam degree of freedom (see below)
      sV <- diag(sireSigma)[i] #sire variance
      dV <- diag(damSigma)[i] #dam variance
      rV <- (diag(indSigma)+scale^2)[i] #residual variance
      ## 
      a.s <- D.s^-1 - sum(D.s^-1)/S
      b.s <- D.s^-1*H.s^-1 - sum(D.s^-1 * H.s^-1)/S
      c.1 <- sum( ((D.s-2)/(sum(D.s)-S)) * tapply((I.sd^-1 - rep(H.s^-1, times=D.s) )^2/rep(D.s,times=D.s) , rep(sort(levels(dat.tren$si)), times=D.s), sum))
      c.2 <- sum( ((D.s-1)/(sum(D.s)-S)) * (H.s^-1 - sum( ((D.s-1)/(sum(D.s)-S))*H.s^-1)  )^2 )
      siDFappREML <- (S-1)*( 1 + S*(S-2)/(S-1) * (sum( (a.s*dV + b.s*rV)^2 ))/((sum( sV + dV*D.s^-1 + rV*D.s^-1*H.s^-1  ))^2)  )^-1
      daDFappREML <- (sum(D.s)-S) * (1 + (rV^2*(c.1+c.2))/(dV + rV*sum( ( (D.s-1)*H.s^-1 )/( sum(D.s)-S ) )^-1)^2)^-1
      siDFappREML.vec <- c(siDFappREML.vec, siDFappREML)
      daDFappREML.vec <- c(daDFappREML.vec, daDFappREML)
### ANOVAuw (with&without unweighted sums of squares)
      if(ANOVApartitioning) {#(anova+++)
        bary.sd <- tapply(dat.tren$y, dat.tren$si_da, mean, na.rm=TRUE)[order(levels(dat.tren$si_da))] ##
        bary.sd <- bary.sd[!is.na(bary.sd)]
        bary.s <- tapply(dat.tren$y, dat.tren$si, mean,na.rm=TRUE)[order(levels(dat.tren$si))] ##
        bary.s <- bary.s[!is.na(bary.s)]
        bary <- mean(dat.tren$y)
        y <- dat[order(dat.tren$sidaid),]$y
        ##
        SS1 <-sum(I.sd*rep((bary.s-bary)^2, times=D.s))
        S1 <- SS1/(S-1)
        SS1U <- w2u*sum((bary.s-bary)^2)
        S1U <- SS1U/(S-1)
        ##
        SS2 <- sum(I.sd*(bary.sd-rep(bary.s, times=D.s))^2)
        S2 <- SS2/daDF
        SS2U <- w3u*sum((bary.sd - rep(bary.s, times=D.s))^2)##
        S2U <- SS2U/daDF
        ##
        SS3 <- sum((y -rep(bary.sd, times=I.sd))^2)
        S3 <- SS3/idDF
        ##
        reV.ANOVA <- S3
        daV.ANOVAuw <- (S2U - S3)/w3u
        daV.ANOVA <- (S2-S3)/mean(I.sd)
        siV.ANOVAuw <- (S1U - (reV.ANOVA + w1u*daV.ANOVAuw))/w2u
        siV.ANOVA <- (S1 - (reV.ANOVA + mean(I.sd)*daV.ANOVA))/(mean(I.sd)*mean(D.s))
        ##
        siDFappANOVA <- (S-1)*( 1 + S*(S-2)/(S-1) * (sum( (a.s*daV.ANOVA + b.s*reV.ANOVA)^2 ))/((sum( siV.ANOVA + daV.ANOVA*D.s^-1 + reV.ANOVA*D.s^-1*H.s^-1  ))^2)  )^-1
        siDFappANOVAuw <- (S-1)*( 1 + S*(S-2)/(S-1) * (sum( (a.s*daV.ANOVAuw + b.s*reV.ANOVA)^2 ))/((sum( siV.ANOVAuw + daV.ANOVAuw*D.s^-1 + reV.ANOVA*D.s^-1*H.s^-1  ))^2)  )^-1
        daDFappANOVA <- (sum(D.s)-S) * (1 + (rV^2*(c.1+c.2))/(daV.ANOVA + rV*sum( ( (D.s-1)*H.s^-1 )/( sum(D.s)-S ) )^-1)^2)^-1
        daDFappANOVAuw <- (sum(D.s)-S) * (1 + (rV^2*(c.1+c.2))/(daV.ANOVAuw + rV*sum( ( (D.s-1)*H.s^-1 )/( sum(D.s)-S ) )^-1)^2)^-1
        ## adding the different estimates to a vector for exporting
        reV.ANOVA.vec <- c(reV.ANOVA.vec, reV.ANOVA)
        daV.ANOVA.vec <- c(daV.ANOVA.vec, daV.ANOVA)
        daV.ANOVAuw.vec <- c(daV.ANOVAuw.vec, daV.ANOVAuw)
        siV.ANOVA.vec <- c(siV.ANOVA.vec, siV.ANOVA)
        siV.ANOVAuw.vec <- c(siV.ANOVAuw.vec, siV.ANOVAuw)
        ##
        siDFappANOVA.vec <- c(siDFappANOVA.vec, siDFappANOVA)
        siDFappANOVAuw.vec <- c(siDFappANOVAuw.vec, siDFappANOVAuw)
        daDFappANOVA.vec <- c(daDFappANOVA.vec, daDFappANOVA)
        daDFappANOVAuw.vec <- c(daDFappANOVAuw.vec, daDFappANOVAuw)
      } #(anova___)
    } #___end of environment-character loop
  }#(unbalanced___)
#################################
  ##output
  if(unbalanced){
    unbal.Est <- list(siDF.vec=siDF.vec, daDF.vec=daDF.vec, idDF.vec=idDF.vec,     siDFappREML.vec=siDFappREML.vec, daDFappREML.vec=daDFappREML.vec,     siDFappANOVA.vec=siDFappANOVA.vec, daDFappANOVA.vec=daDFappANOVA.vec,     siDFappANOVAuw.vec=siDFappANOVAuw.vec, daDFappANOVAuw.vec=daDFappANOVAuw.vec,      w1u.vec=w1u.vec, w2u.vec=w2u.vec, w3u.vec=w3u.vec) #degrees of freedom and weights
  }else{
    unbal.Est <- list()
  }
  ##
    if(ANOVApartitioning & REMLpartitioning & chN*enN==1){
      secondpartitioning.reml <- list(ANOVAuw=list(para=new("para", rbS=matrix(0), siS=matrix(siV.ANOVAuw), daS=matrix(daV.ANOVAuw), idS=matrix(reV.ANOVA), phS=matrix(), error=0, fixe = array(0, dim=c(enN, chN, fbN)))),
                                        ANOVA=new("para", rbS=matrix(0), siS=matrix(siV.ANOVA), daS=matrix(daV.ANOVA), idS=matrix(reV.ANOVA), phS=matrix(), error=0, fixe = array(0, dim=c(enN, chN, fbN)))
                                      )
    }else{
      secondpartitioning.reml <- list()
    }
#################################
### OUTPUT
  ## orig
  origEst <- new("orig", hist=c(paraDATA@orig@hist, "est"), warn=c(paraDATA@orig@warn, ""), time=c(paraDATA@orig@time, date()), part=paraDATA@orig@part)
  ## supl
  suplEst <- new("supl", paraDATA@supl)
  ## para & spec
  if(part=="ANOVA" & chN*enN==1){
    paraEst <- new("para", rbS=matirx(0), siS=matrix(siV.ANOVA), daS=matrix(daV.ANOVA), idS=matrix(reV.ANOVA), phS=matrix(), error=0, fixe = rep(0, times=chN*enN))
    specEst <- new("spec", additional.partitioning=list(), unbalanced=unbal.Est, modelsummary=list(), secondcontrast=list())
  }
  if(part=="ANOVAuw" & chN*enN==1){
    paraEst <- new("para", rbS=matirx(0), siS=matrix(sV.ANOVAuw), daS=matrix(daV.ANOVAuw), idS=matrix(reV.ANOVAuw), phS=matrix(), error=0, fixe = rep(0, times=chN*enN))
    specEst <- new("spec", additional.partitioning=list(), unbalanced=unbal.Est, modelsummary=list(), secondcontrast=list())
  }
  if(part=="REML"|part=="all"){
    paraEst <- reml.para
    #specEst <- new("spec", list(unbal=unbal, modelsummary=modelsummary.reml, secondcontrast=secondcontrast.reml, secondpartitioning=secondpartitioning.reml))
    specEst <- new("spec", additional.partitioning=secondpartitioning.reml, unbalanced=unbal.Est, modelsummary=modelsummary.reml, secondcontrast=secondcontrast.reml)
  }
  ## DATA
  DATAEst <- new("DATA")
  ## paraDATA
  est <- new("paraDATA", orig=origEst, supl=suplEst, para=paraEst, DATA=DATAEst, spec=specEst)
  ##
  if (file){
    dir.create(path=path, showWarnings = FALSE, recursive = TRUE)
    save(est, file=paste(path, "est.rda", sep=""))
  }
  ##
  est
}
