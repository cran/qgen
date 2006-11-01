sim <- function(paraDATA, file=FALSE){
  ## to d:o
  ## remove design
  ## all that goes into class "supl" (miss -> realized variance-covariance matrix, fixe)
##################################################
  ## helping:
  ## give also the random shocks
  randomshock <- FALSE
  ##
##################################################
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
  sS.pop = paraDATA@para@siS
  dS.pop = paraDATA@para@daS
  iS.pop = paraDATA@para@idS
  error = paraDATA@para@error
  fix = paraDATA@para@fixe
############################################################
#                                                      INTRO
# loading libraries and functions
#
  zeroVector <- rep(0, times=chN*enN)                     # zero vector for means of the varcov matrix       
  zeroMatrix <- matrix(0, ncol=chN*enN, nrow=chN*enN)     # Matrix with all zeros of dimension of the fixed effects
#
############################################################
#                                                    1. PART
# generates the "Designmatrix"
#
  en <- factor(rep(dimnames(fix)[[1]], each=fbN*rbN*siN*daN*idN*design[7]*design[8], times=chN))
  ch <- factor(rep(dimnames(fix)[[2]], each=fbN*enN*rbN*siN*daN*idN*design[7]*design[8]))
  fb <- factor(rep(dimnames(fix)[[3]], each=fbN*siN*daN*idN*design[7]*design[8], times=chN*enN))
  rb <- factor(rep(1:rbN, each=fbN*siN*daN*idN*design[7]*design[8], times=chN*enN))##work space!!
  si <- factor(rep(rep(1:siN, each=((enN*rbN*daN*idN*design[7]*design[8])/(enN*rbN)),times =fbN*rbN), times=chN*enN))
  da <- factor(rep(rep(1:(siN*daN), each=(idN*design[7]*design[8]), times = fbN*rbN), times=chN*enN))
  id <- factor(rep(1:(enN*fbN*rbN*siN*daN*idN), each=design[7]*design[8], times=chN))
  ##
  fm <- data.frame(ch, en, fb, rb, si, da, id)
  fm.sorted <- fm[order(fm$en),] #damit zuerst alle en1..
##########################################################
###                                                  2. PART
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
############################################################
  cc <- cbind(fm.sorted, dc)
############################################################
#                                                    3. PART
# calculation of y
#
  y <- as.vector(fix[cbind(unclass(cc$en),unclass(cc$ch),unclass(cc$fb))] + cc$randomblockcov + cc$sirecov + cc$damcov + cc$indcov + cc$errorvar) #as.vector is IMPORTANT
  data <- cbind(cc,y)
  ## removes the random shocks from "data" if not explicitely wished
  if(!randomshock){
    data$randomblockcov <- NULL
    data$sirecov <- NULL
    data$damcov <- NULL
    data$indcov <- NULL
    data$errorvar <- NULL
  }
############################################################
#                                                    4. PART
#     UNBALANCEDNESS
  exclude <- round(miss*rbN*siN*daN*idN*design[7]*design[8],0)
  fail <- NULL
  for(i in 1:length(exclude)){
    fail <- c(fail, sample((((i-1)*rbN*siN*daN*idN*design[7]*design[8])+1):(i*rbN*siN*daN*idN*design[7]*design[8]), exclude[i], replace=FALSE))
  }
  if(length(fail)>1){
    data <- data[-fail,]
  }
  ### realized Variance
  # not yet implemented for MISSING values (idea take data with shocks; minus fail; take
##     en_ch = factor(paste(en,ch, sep="_"))
##   if (rbN>1){
##     randomblockshock.realized <- table(data$si, data$en_ch)!=0### * randomblockshock
##     tS.rea <- cov(randomblockshock.realized, use="pairwise.complete.obs")
##   } else {
##     tS.rea <- zeroMatrix
##   }
##   if (siN>1){
##     sipresent <- table(data$si, data$en_ch)!=0 # a table with FALSE/TRUE -> numeric:0/1 for presence absence
##     sireshock.realized <- sireshock * sipresent
##     sS.rea <- cov(sireshock.realized, use="pairwise.complete.obs")
##   } else {
##     sS.rea <- zeroMatrix
##   }
##   if (daN>1){
##     dapresent <- table(data$da, data$en_ch)!=0 # a table with FALSE/TRUE -> numeric:0/1 for presence absence
##     damshock.realized <- damshock * dapresent
##     dS.rea <- cov(damshock.realized, use="pairwise.complete.obs")
##   } else {
##     dS.rea <- zeroMatrix
##   }
##   if (rbN*siN*daN > 1){
##     trueid <- factor(rep(1:(fbN*rbN*siN*daN*idN), each=design[7]*design[8], times=chN*enN))
##     idpresent <- table(data$trueid, data$en_ch)!=0
##     indshock.realized <- indshock * idpresent
##     iS.rea <- cov(indshock.realized, use="pairwise.complete.obs")
##   } else {
##     iS.rea <- zeroMatrix
##   }
#   entr.means <- tapply(data$y, data$en_ch, mean)
#   entr.means <- matrix(en_ch.means, nrow=1, dimnames=list(NULL, names(entr.means)))
#   tr.means <- tapply(data$y, data$ch, mean)
#   tr.means <- matrix(t.means, nrow=1, dimnames=list(NULL, names(tr.means)))
#   en.means <- tapply(data$y, data$en, mean)
#   en.means <- matrix(en.means, nrow=1, dimnames=list(NULL, names(en.means)))
############################################################
  ##                                                    5. PART
  ##     defines the structure of the output                    
  ##
  origSim <- new("orig", hist=c(paraDATA@orig@hist, "sim"), warn=c(paraDATA@orig@warn, ""), time=c(paraDATA@orig@time, date()), part=paraDATA@orig@part)
  suplSim <- new("supl", paraDATA@supl)
  #paraSim <- new("para", rbS=tS.rea, siS=sS.rea, daS=dS.rea, idS=iS.rea, phS=matrix(), error=error, fixe=fix)                   #fix needs to be programed!!!
  paraSim <- new("para")
  DATASim <- new("DATA", dat=data)
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
############################################################
#                                                        END
