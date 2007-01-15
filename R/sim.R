#paraDATA <- the(enN=2, chN=3, fbN=1, siN=3, daN=3, idN=2)# miss=array(c(0.1,0.9,0.3,0.7), dim=c(5,4,1)))




sim <- function(paraDATA, file=FALSE){
  ## to d:o
  ## remove design
  ## all that goes into class "supl" (miss -> realized variance-covariance matrix, fixe)
##################################################
  ## helping:
  ##
#############################################################################
###                                                                   READING
###
  chN = paraDATA@supl@chN
  enN = paraDATA@supl@enN
  fbN = paraDATA@supl@fbN
  rbN = paraDATA@supl@rbN
  siN = paraDATA@supl@siN
  daN = paraDATA@supl@daN
  idN = paraDATA@supl@idN
  design=c(0,0,0,0,0,0,1,1)# to be replaced in a later version (would allow further nested random variables)
  miss = paraDATA@supl@miss
  ##
  rbS.pop = paraDATA@para@rbS
  sS.pop = paraDATA@para@siS  # also as source for the dimname of covariance matrices
  dS.pop = paraDATA@para@daS
  iS.pop = paraDATA@para@idS
  error = paraDATA@para@error
  fix = paraDATA@para@fixe    # source for the names of the fixed effects
#############################################################################
###                                                                     INTRO
###
  ## names
  S.names <- dimnames(sS.pop)[[1]]     # the names of rows and columns in variance covariance matrix
  ##
  zeroVector <- rep(0, times=chN*enN)                     # zero vector for means of the varcov matrix       
  zeroMatrix <- matrix(0, ncol=chN*enN, nrow=chN*enN, dimnames=list(S.names,S.names))   # Matrix with all zeros of dimension of the fixed effects
#
#############################################################################
###                                                             DESIGN MATRIX
###
  en <- factor(rep(dimnames(fix)[[1]], each=fbN*rbN*siN*daN*idN*design[7]*design[8], times=chN))
  ch <- factor(rep(dimnames(fix)[[2]], each=fbN*enN*rbN*siN*daN*idN*design[7]*design[8]))
  fb <- factor(rep(dimnames(fix)[[3]], each=fbN*siN*daN*idN*design[7]*design[8], times=chN*enN))
  rb <- factor(rep(1:rbN, each=fbN*siN*daN*idN*design[7]*design[8], times=chN*enN))##work space!!
  si <- factor(rep(rep(1:siN, each=((enN*rbN*daN*idN*design[7]*design[8])/(enN*rbN)),times =fbN*rbN), times=chN*enN))
  da <- factor(rep(rep(1:(siN*daN), each=(idN*design[7]*design[8]), times = fbN*rbN), times=chN*enN))
  id <- factor(rep(1:(enN*fbN*rbN*siN*daN*idN), each=design[7]*design[8], times=chN))
  thid <- factor(rep(1:(fbN*rbN*siN*daN*idN), each=design[7]*design[8], times=chN*enN))# theoretical individual; is used to calculate the realized Sigma for level individual! 
  ##
  fm <- data.frame(ch, en, fb, rb, si, da, id, thid)
  fm.sorted <- fm[order(fm$en),] #damit zuerst alle en1..
#############################################################################
###                                                              RANDOM SHOCK
###
###     generates random shocks (mean 0, cov-matrix as demanded) for all random effects
###
  ##___randomblockCOV
  if(rbN>1){
    randomblockshock <- mvrnorm(rbN, mu=zeroVector, Sigma=rbS.pop)
    randomblockcov <- rep(c(randomblockshock), each=fbN*siN*daN*idN*design[7]*design[8])
  }else{randomblockshock <- NULL; randomblockcov <- 0}
  ##___sireCOV
  sireshock <- mvrnorm(siN, mu=zeroVector, Sigma=sS.pop)
  sirecov <- rep(c(sireshock[, rep(1:(chN*enN), each=fbN*rbN)]), each=(daN*idN*design[7]*design[8]))
  ##___damCOV
  damshock <- mvrnorm(siN*daN, mu=zeroVector, Sigma=dS.pop)
  damcov <- rep(c(damshock[, rep(1:(chN*enN), each=fbN*rbN)]), each=(idN*design[7]*design[8]))
  ##___indCOV  for every individual within all trait x environmnet combinations
  indshock <- mvrnorm(fbN*rbN*siN*daN*idN*design[7]*design[8], mu = zeroVector, Sigma=iS.pop)
  indcov <- c(indshock)
  ##___errorCOV  for every individual within all trait x environmnet combinations
  errorshock <- rnorm(fbN*rbN*siN*daN*idN*design[7]*design[8] ,mean = 0, sd=sqrt(error))
  errorvar <- c(errorshock)
  ##
  dc <- data.frame(randomblockcov, sirecov, damcov, indcov, errorvar)
  ##
  cc <- cbind(fm.sorted, dc)
#############################################################################
###                                                          CALCULATION OF Y
###
  y <- as.vector(fix[cbind(unclass(cc$en),unclass(cc$ch),unclass(cc$fb))] + cc$randomblockcov + cc$sirecov + cc$damcov + cc$indcov + cc$errorvar) #as.vector is IMPORTANT
  data <- cbind(cc,y)
#############################################################################
###                                                            UNBALANCEDNESS
###
### miss is a vector length 3 (fb, en, ch) (see paraDATA@supl@miss)
  exclude <- round(miss*rbN*siN*daN*idN*design[7]*design[8],0) # how many measurements are missing? 
  fail <- NULL
  for(i in 1:length(exclude)){
    fail <- c(fail, sample((((i-1)*rbN*siN*daN*idN*design[7]*design[8])+1):(i*rbN*siN*daN*idN*design[7]*design[8]), exclude[i], replace=FALSE))
  }
  if(length(fail)>1){
    data <- data[-fail,]
  }
#############################################################################
###                                                         REALIZED VARIANCE
### takes the ..shocks and calculates the variance-covariance matrix
### if some shocks were not used because of missing values they are excluded
  data$ench = factor(paste(data$en,data$ch, sep="_"))
  ##
  if (rbN>1){
    rb.shock.realized <- table(data$si, data$fbench)!=0
    rbS.rea <- cov(rb.shock.realized, use="pairwise.complete.obs")
  } else {
    rbS.rea <- zeroMatrix
  }
  ##
  if (siN>1){
    si.present <- table(data$si, data$ench)!=0  # a table with FALSE/TRUE -> numeric:0/1 for presence absence
    si.shock.realized <- sireshock * si.present  # removes the realized shockes that are not present at all because all individuals are missing
    siS.rea <- cov(si.shock.realized, use="pairwise.complete.obs")
  } else {
    siS.rea <- zeroMatrix
  }
  ##
  if (daN>1){
    da.present <- table(data$da, data$ench)!=0 # a table with FALSE/TRUE -> numeric:0/1 for presence absence
    da.shock.realized <- damshock * da.present
    daS.rea <- cov(da.shock.realized, use="pairwise.complete.obs")
  } else {
    daS.rea <- zeroMatrix
  }
  ## individuals
  ## because every environment has its own individuals we have the factor
  ##    "thid" which has a level for individuals but not for those which
  ##    are in a new environment (see its definition above)
  if (rbN*siN*daN > 1){
    id.present <- table(data$thid, data$ench)
    id.shock.realized <- indshock * id.present
    idS.rea <- cov(id.shock.realized, use="pairwise.complete.obs")
  } else {
    idS.rea <- zeroMatrix
  }
#   entr.means <- tapply(data$y, data$en_ch, mean)
#   entr.means <- matrix(en_ch.means, nrow=1, dimnames=list(NULL, names(entr.means)))
#   tr.means <- tapply(data$y, data$ch, mean)
#   tr.means <- matrix(t.means, nrow=1, dimnames=list(NULL, names(tr.means)))
#   en.means <- tapply(data$y, data$en, mean)
#   en.means <- matrix(en.means, nrow=1, dimnames=list(NULL, names(en.means)))
#############################################################################
###                                                                    OUTPUT
### 
### defines the structure of the output
  origSim <- new("orig", hist=c(paraDATA@orig@hist, "sim"), warn=c(paraDATA@orig@warn, ""), time=c(paraDATA@orig@time, date()), part=paraDATA@orig@part)
  suplSim <- new("supl", paraDATA@supl)
  paraSim <- new("para", rbS=rbS.rea, siS=siS.rea, daS=daS.rea, idS=idS.rea, phS=matrix(), error=error, fixe=fix)                   #fix needs to be programed!!!
  DATASim <- new("DATA", dat=data.frame(ch=data$ch, en=data$en, fb=data$fb, rb=data$rb, si=data$si, da=data$da, id=data$id, y=data$y))
  specSim <- new("spec")
  ##
  sim <- new("paraDATA", orig=origSim, supl=suplSim, para=paraSim, DATA=DATASim, spec=specSim)
  ##
  if (file){
    dir.create(path=path, showWarnings = FALSE, recursive = TRUE)
    save(sim, file=paste(path, "sim.rda", sep=""))
  }
  sim
}
#############################################################################
###                                                                       END
# known: fix is just what it got from within the paraDATA on average it will not have changed but only on average...
