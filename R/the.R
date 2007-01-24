the <- function(partitioning = "all",
                chN = 1, # character
                enN = 1, # environment
                fbN = 1, # fixeblock
                rbN = 1, # randomblock
                siN = 100, # sire
                daN = 6, # dam
                idN = 3, # individuals
                randomblockCor = matrix(0.5, chN*enN, chN*enN),
                randomblockVar = rep(100, chN*enN),
                additiveCor = matrix(0.5, chN*enN, chN*enN),
                additiveVar = rep(100, chN*enN),
                dominanceCor = matrix(0.5, chN*enN, chN*enN),
                dominanceVar = rep(100, chN*enN),
                maternalCor = matrix(0.5, chN*enN, chN*enN),
                maternalVar = rep(100, chN*enN),
                environmentalCor = matrix(0, chN*enN, chN*enN), # it MUST be zero if only one measurement per individual in every environment!
                environmentalVar = rep(100, chN*enN),
                ch.names = paste("ch",1:chN,sep=""),
                en.names = paste("en",1:enN,sep=""),
                fb.names = paste("fb",1:fbN,sep=""),
                fixe = array(0, dim=c(enN, chN, fbN)),
                miss = array(0, dim=c(enN, chN, fbN)),
                file=TRUE,
                path="~/qgen/")
### the arguments:
#                 partitioning = "all"
#                 chN = 1, # character
#                 enN = 1, # environment
#                 fbN = 1, # fixeblock
#                 rbN = 1, # randomblock
#                 siN = 100, # sire
#                 daN = 6, # dam
#                 idN = 3, # individuals
#                 randomblockCor = matrix(0.5, chN*enN, chN*enN);
#                 randomblockVar = rep(0, chN*enN);
#                 additiveCor = matrix(0.66, chN*enN, chN*enN);
#                 additiveVar = rep(100, chN*enN);
#                 dominanceCor = matrix(0.66, chN*enN, chN*enN);
#                 dominanceVar = rep(100, chN*enN);
#                 maternalCor = matrix(0.66, chN*enN, chN*enN);
#                 maternalVar = rep(100, chN*enN);
#                 environmentalCor = matrix(0, chN*enN, chN*enN); # it MUST be zero if only one measurement per individual in every environment!
#                 environmentalVar = rep(100, chN*enN);
#                 ch.names = paste("ch",1:chN,sep="");
#                 en.names = paste("en",1:enN,sep="");
#                 fb.names = paste("fb",1:fbN,sep="");
#                 fixe = array(0, dim=c(enN, chN, fbN));
#                 miss = array(0, dim=c(enN, chN, fbN));
#                 file=TRUE
{
### missing (transformed into an array if only a scalar is given)
  miss <- array(miss, dim=c(enN, chN, fbN))
### SIGMA construction
### building a variance-covariance matrix from correalation matrix and variance vector
  SigmaConstructor <- function(var,corr){
    n <- dim(corr)[1]
    corr[lower.tri(corr,diag=TRUE) & upper.tri(corr, diag = TRUE)] <- 1 #sets diagonal elements to exactly 1
    Sigma <- corr * matrix(sqrt(var), nrow=n, ncol=n) * matrix(sqrt(var), nrow=n, ncol=n, byrow=TRUE)
    Sigma
  }
### CAUSUAL variance-covariance matrices
  if(rbN>1){randomblockSigma <- SigmaConstructor(var=randomblockVar, cor=randomblockCor)}else{randomblockSigma <- NULL}
  additiveSigma <- SigmaConstructor(var=additiveVar, cor=additiveCor)
  dominanceSigma <- SigmaConstructor(var=dominanceVar, cor=dominanceCor)
  maternalSigma <- SigmaConstructor(var=maternalVar, cor=maternalCor)
  environmentalSigma <- SigmaConstructor(var=environmentalVar, cor=environmentalCor)
### OBSERVABLE variance-covariance matrices
  ## from underlying causal components
  zeroVector <- rep(0, times=chN*enN)                     # zero vector for means of the varcov matrix       
  zeroMatrix <- matrix(0, ncol=chN*enN, nrow=chN*enN)     # Matrix with all zeros of dimension of the fixed effects
  covarianceDimNames <- list(paste(rep(en.names,each=chN),rep(ch.names,times=enN),sep="_"), paste(rep(en.names,each=chN),rep(ch.names,times=enN),sep="_"))
  ##
  if (rbN>1){
    rbS <- matrix(randomblockSigma, dimnames=covarianceDimNames, nrow=chN*enN, ncol=chN*enN)
  }else{rbS <- matrix(zeroMatrix, dimnames=covarianceDimNames, nrow=chN*enN, ncol=chN*enN)}
  siS <- matrix(1/4*additiveSigma, dimnames=covarianceDimNames, nrow=chN*enN, ncol=chN*enN)
  daS <- matrix(1/4*additiveSigma + 1/4*dominanceSigma + maternalSigma,              dimnames=covarianceDimNames, nrow=chN*enN, ncol=chN*enN)    
  iS <- matrix(1/2*additiveSigma + 3/4*dominanceSigma + environmentalSigma,         dimnames=covarianceDimNames, nrow=chN*enN, ncol=chN*enN)     
  idS <- iS * kronecker(diag(1,nrow=enN, ncol=enN), matrix(1, nrow=chN,ncol=chN)) # makes the individual Sigma having the OBSERVABLE form 
  phS <- matrix(additiveSigma + dominanceSigma + maternalSigma + environmentalSigma, dimnames=covarianceDimNames, nrow=chN*enN, ncol=chN*enN) # allows the CORRECT calculation of heritability of plasticity but can NOT be observed
  error <- 0
### para_fixe
  dimnames(fixe) <- list(en.names, ch.names, fb.names)
  dimnames(miss) <- list(en.names, ch.names, fb.names)
### OUTPUT
  origTHE <- new("orig", hist="the", warn="", time=date(), part=partitioning)
  suplTHE <- new("supl", chN=as.integer(chN), enN=as.integer(enN), fbN=as.integer(fbN), rbN=as.integer(rbN), siN=as.integer(siN), daN=as.integer(daN), idN=as.integer(idN), miss=miss)
  fixe.out <- fixe
  paraTHE <- new("para", rbS=rbS, siS=siS, daS=daS, idS=idS, phS=phS, error=0, fixe=fixe.out)
  DATATHE <- new("DATA")
  specTHE <- new("spec")
  ##
  the <- new("paraDATA", orig=origTHE, supl=suplTHE, para=paraTHE, DATA=DATATHE, spec=specTHE)
  ##
  if (file){
    dir.create(path=path, showWarnings = FALSE, recursive = TRUE)
    save(the, file=paste(path, "the.rda", sep=""))
  }
  the
}
