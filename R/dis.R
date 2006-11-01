dis <- function(path="~/Desktop/myproject.qgen/", alpha=0.05){
### THE T LEVEL
  if(file.exists(paste(path,"statT.rda",sep=""))){
    input.name <- load(file=paste(path,"statT.rda",sep=""))
    eval(parse(text=paste("theta", " <- ", input.name)))
    rm(list=c(input.name))
    rm(input.name) 
  }
### THE S LEVEL
  if(file.exists(paste(path,"statS.rda",sep=""))){
    input.name <- load(file=paste(path,"statS.rda",sep=""))
    eval(parse(text=paste("t.hat", " <- ", input.name)))
    rm(list=c(input.name))
    rm(input.name)
    ##
    S.number <- dim(t.hat)[1]
    ##
    t.hat.S <- array(t.hat, dim=c(S.number), dimnames=list(substring(dimnames(t.hat)[[1]], first=3, last=7)))
  }
### THE R LEVEL
  if(file.exists(paste(path,"statR.rda",sep=""))){
    input.name <- load(file=paste(path,"statR.rda",sep=""))
    eval(parse(text=paste("t.star", " <- ", input.name)))
    rm(list=c(input.name))
    rm(input.name)
    ##
    R.number <- table(substring(dimnames(t.star)[[1]], first=3, last=7))[1]
    ##
    t.hat.SR <- array(t.hat, dim=c(S.number,R.number), dimnames=list(paste("S",leading(1:S.number),sep=""),paste("R",leading(1:R.number),sep="")))
    ##
    t.star.SR <- tapply(t.star, list(substring(dimnames(t.star)[[1]], first=3, last=7),substring(dimnames(t.star)[[1]], first=8, last=12)), function(x){x})
    ## sorts t.star.SR
    t.star.SR.sorted <- t.star.SR
    for (i in seq(1:dim(t.star.SR.sorted)[1])){
      t.star.SR.sorted[i,] <- sort(t.star.SR.sorted[i,], na.last=TRUE)
    }
    ##
    bias.star.SR <- t.star.SR - t.hat.SR
    ##
    v.star.SR <- array(tapply(t.star, substring(dimnames(t.star)[[1]], first=3, last=7), var, na.rm=TRUE), dim=c(S.number, R.number), dimnames=list(paste("S",leading(1:S.number),sep=""),paste("R",leading(1:R.number),sep="")))
    v.star.S <- tapply(t.star, substring(dimnames(t.star)[[1]], first=3, last=7), var, na.rm=TRUE)
    ## CI Percentile
    perc.L.star.S <- interpolation(R=R.number, alpha=alpha, sort.mat=t.star.SR.sorted)
    perc.U.star.S <- interpolation(R=R.number, alpha=1-alpha, sort.mat=t.star.SR.sorted)
    ## CI basic
    basic.L.star.S <- 2*t.hat.S - interpolation(R=R.number, alpha=1-alpha, sort.mat=t.star.SR.sorted)
    basic.U.star.S <- 2*t.hat.S - interpolation(R=R.number, alpha=alpha, sort.mat=t.star.SR.sorted)
    ##  BCa
#    m3 <- colSums((t.star.S - matrix(colMeans(t.star.S, na.rm=TRUE), nrow=dim(t.star)[1], ncol=dim(t.star)[2],byrow=TRUE))^3, na.rm=TRUE) / dim(t.star)[1] 
#    s3 <- sqrt(colVars(t.star, na.rm=TRUE))^3
#    a = (m3/s3)/6
#    w = qnorm(colSums(ifelse(t.star <= rep(t.hat, R.number), 1, 0),na.rm=TRUE)/(R.number+1))
#     z.tilde.alpha.low = w + qnorm(alpha)
#     alpha.tilde.low = pnorm(w + z.tilde.alpha.low/(1 - a*z.tilde.alpha.low))
#     z.tilde.alpha.high= w + qnorm(1-alpha)
#     alpha.tilde.high = pnorm(w + z.tilde.alpha.high/(1 - a*z.tilde.alpha.high))
#     BCa.L.star <- tapply(t.star, substring(dimnames(t.star)[[1]], first=3, last=7), function(x){interpolation(R=length(x), alpha=alpha.tilde.low, sort.mat=as.matrix(sort(x)))})
#     BCa.U.star <- tapply(t.star, substring(dimnames(t.star)[[1]], first=3, last=7), function(x){interpolation(R=length(x), alpha=alpha.tilde.high, sort.mat=as.matrix(sort(x)))})
   }
### THE Q LEVEL
  if(file.exists(paste(path,"statQ.rda",sep=""))){
    input.name <- load(file=paste(path,"statQ.rda",sep=""))
    eval(parse(text=paste("t.starstar", " <- ", input.name)))
    rm(list=c(input.name))
    rm(input.name)
    ##
    v.starstar.SR <- tapply(t.starstar, list(substring(dimnames(t.starstar)[[1]], first=3, last=7),substring(dimnames(t.starstar)[[1]], first=8, last=12)),var)
    ## sorts z.star.SR: each S individually
    z.star.SR.sorted <- (t.star.SR - t.hat.SR) / sqrt(v.starstar.SR)    
    for (i in seq(1:dim(z.star.SR.sorted)[1])){
      z.star.SR.sorted[i,] <- sort(z.star.SR.sorted[i,], na.last=TRUE)
    }
    dimnames(z.star.SR.sorted)[[2]] <- rep("Rsorted", times=dim(z.star.SR.sorted)[2])
    ## CI studentized
    stud.L.star.S <- t.hat.S - sqrt(v.star.S) * interpolation(R=R.number, alpha=1-alpha, sort.mat=z.star.SR.sorted)
    stud.U.star.S <- t.hat.S - sqrt(v.star.S) * interpolation(R=R.number, alpha=alpha, sort.mat=z.star.SR.sorted)
  }
  ci <- function(...){
    print(paste("percentile :      ", perc.L.star.S," - ", perc.U.star.S))
    print(paste("basic      :      ", basic.L.star.S," - ", basic.U.star.S))
    #print(paste("studentized:      ", stud.L.star.S," - ", stud.U.star.S))
  }
  ci()
}



#   ## transformation="none"; alpha=0.05; stat="HER"; dPLOT=TRUE; sPLOT=TRUE; est.ci=TRUE; presentation=TRUE
#   ## output: "tex", "sPLOTdata", "sPLOTbw" "sPLOTcolor" "dPLOT" "CIlength"
#   ## statistic=Stat2e
#   ## load multi
#   load(paste(path,folder,"/","multi.rda", sep=""))
#   ## reading general statistic properties
  
#   statlimes <- attributes(statistic)$limes
# ### POP level
#   missing <- multi$pop$supl$mis



#                                         #design <- multi$pop$supl$design
#   S <- 0; R <- 0; Q <- 0; N <- NULL

#   ##  theta
#   theta <- tnf(statistic(multi$pop, alpha=alpha, frommethod=frommethod)$spec$HER)
#   ##
#   S <- length(multi$call)
#   N <- NULL; studL <- NULL; studU <- NULL; percL <- NULL;  percU <- NULL;  BCaL <- NULL;  BCaU <- NULL;  basicL <- NULL;  basicU <- NULL;  unwL <- NULL;  unwU <- NULL;  t.hat <- NULL;  senL <- NULL;  senU <- NULL;  t.rea <- NULL
# ### loop S
#   for (s in 1:S){ # S loop
#     if (!is.null(multi$call[[s]]$est)){ # S cond
#       N.s <- NULL
#       t.hatL.s <- NULL
#       t.hatU.s <- NULL
#       t.rea.s <- tnf(statistic(multi$call[[s]]$rea, alpha=alpha, frommethod=frommethod)$spec$HER)
#       t.hat.s <- tnf(statistic(multi$call[[s]]$est, alpha=alpha, frommethod=frommethod)$spec$HER)
#       senU <- rbind(senU, tnf(statistic(multi$call[[s]]$est, alpha=alpha, frommethod=frommethod)$spec$Upper))
#       senL <- rbind(senL, tnf(statistic(multi$call[[s]]$est, alpha=alpha, frommethod=frommethod)$spec$Lower))
#       ##
#       R <- length(multi$call[[s]]$boot)
# #### cond R
# 		if (R>0){ ####cond R     
# 		    t.star <- NULL; v.star <- NULL;  z.star <- NULL
#                     for (i in 1:R){#### loop R
# 			statREA.i <- NULL; t.star.i <- NULL; t.starstar <- NULL; v.i <- NULL; z.star.i <- NULL
#                         if (stat=="booot") t.star.i <- tnf(booot(multi$call[[s]]$boot[[i]]$star))
# 			if (stat=="HER") t.star.i <- tnf(Stat1(multi$call[[s]]$boot[[i]]$star, alpha=alpha)$meth$HER)
# 			if (stat=="HER_unw") t.star.i <- tnf(Stat1unw(multi$call[[s]]$boot[[i]]$star, alpha=alpha)$meth$HER_unw)
# 			if(stat=="HER_ind") t.star.i <- tnf(Stat2e(multi$call[[s]]$boot[[i]]$star)$meth$HER_ind)
# 			if(stat=="aC")   {aC.i <- as.matrix(StatXt(multi$call[[s]]$boot[[i]]$star)$meth$aC);       t.star.i <- tnf(matrix(aC.i[lower.tri(aC.i)],nrow=1))}
# 			if(stat=="aC.i") {aC.i.i <- as.matrix(Stat2e(multi$call[[s]]$boot[[i]]$star)$meth$sC.i);   t.star.i <- tnf(matrix(aC.i.i[lower.tri(aC.i.i)],nrow=1)) }
# 			if(stat=="aC.ic"){aC.ic.i <- as.matrix(Stat2e(multi$call[[s]]$boot[[i]]$star)$meth$sC.ic); t.star.i <- tnf(matrix(aC.ic.i,nrow=1))}
# 			if(stat=="aC.ih"){aC.ih.i <- as.matrix(Stat2e(multi$call[[s]]$boot[[i]]$star)$meth$sC.ih); t.star.i <- tnf(matrix(aC.ih.i,nrow=1))}
# 			#
# 			t.star <- rbind(t.star, t.star.i)
# 			#
# 			t.starstar.j <- NULL
# 			Q <- length(multi$call[[s]]$boot[[i]]$starstar)
# 			N.s <- c(N.s , Q)
# 			if (Q>0){#### cond Q
# 			    for (j in 1:Q){#### loop Q
#                               	if (stat=="booot") t.starstar.j <- tnf(booot(multi$call[[s]]$boot[[i]]$starstar[[j]]))
# 				if (stat=="HER") t.starstar.j <- tnf(Stat1(multi$call[[s]]$boot[[i]]$starstar[[j]], alpha=alpha)$meth$HER)
# 				if (stat=="HER_unw") t.starstar.j <- tnf(Stat1unw(muli$call[[s]]$boot[[i]]$starstar[[j]], alpha=alpha)$meth$HER_unw)
# 				if (stat=="HER_ind") t.starstar.j <- tnf(Stat2e(multi$call[[s]]$boot[[i]]$starstar[[j]])$meth$HER_ind)
# 				if (stat=="aC")   {aC.j <- as.matrix(StatXt(multi$call[[s]]$boot[[i]]$starstar[[j]])$meth$aC);       t.starstar.j <- tnf(matrix(aC.j[lower.tri(aC.j)],nrow=1))}
# 				if (stat=="aC.i") {aC.i.j <- as.matrix(Stat2e(multi$call[[s]]$boot[[i]]$starstar[[j]])$meth$sC.i);   t.starstar.j <- tnf(matrix(aC.i.j[lower.tri(aC.i.j)],nrow=1))}
#                                 if (stat=="aC.ic"){aC.ic.j <- as.matrix(Stat2e(multi$call[[s]]$boot[[i]]$starstar[[j]])$meth$sC.ic); t.starstar.j <- tnf(matrix(aC.ic.j,nrow=1))}
# 				if (stat=="aC.ih"){aC.ih.j <- as.matrix(Stat2e(multi$call[[s]]$boot[[i]]$starstar[[j]])$meth$sC.ih); t.starstar.j <- tnf(matrix(aC.ih.j,nrow=1))}
# 				#
# 				t.starstar <- rbind(t.starstar, t.starstar.j)
# 			    }#____ loop Q
# 			    v.i <- colVars(t.starstar, na.rm=TRUE)
# 			    v.star <- rbind(v.star, v.i)
# 			    z.star.i <- (t.star.i - t.hat.s)/sqrt(v.i)
# 			    z.star <- rbind(z.star, z.star.i)
# 			} #_____ cond Q
# 		    } #_________________________ loop R
# 		    ## BIAS and VARIANCE
# 		    t.hat.mat <- matrix(t.hat.s, nrow=R, ncol=dim(t.hat.s)[2], byrow=TRUE)
# 		    bias <- t.star - t.hat.mat # if theta is known there is a better "bias estimator" available see below
# 		    v <- colVars(t.star)
# 		    ## Preparation for CI
# 		    t.star.sort <- NULL
# 		    for (i in (1:dim(t.star)[2])){t.star.sort <- cbind(t.star.sort, sort(t.star[,i]))}
#                     ## CI studentized
# 		    if (all(N.s>4)) {
#                       z.star.sort <- NULL
#                       for (i in (1:dim(z.star)[2])){z.star.sort <- cbind(z.star.sort, sort(z.star[,i]))}
#                       studL <- rbind(studL, t.hat.s - sqrt(v) * interpolation(R=R, alpha=1-alpha, sort.mat=z.star.sort))
#                       studU <- rbind(studU, t.hat.s - sqrt(v) * interpolation(R=R, alpha=alpha, sort.mat=z.star.sort))
#                     }else{studL<- rbind(studL, rep(NA,dim(t.star)[2])); studU <- rbind(studU, rep(NA,dim(t.star)[2]))}
#                     ## CI Percentile
#                     percL <- rbind(percL, interpolation(R=R, alpha=alpha, sort.mat=t.star.sort))
#                     percU <- rbind(percU, interpolation(R=R, alpha=1-alpha, sort.mat=t.star.sort))
#                     ##  BCa
#                     m3 <- colSums((t.star - matrix(colMeans(t.star, na.rm=TRUE),nrow=dim(t.star)[1], ncol=dim(t.star)[2],byrow=TRUE))^3, na.rm=TRUE) / dim(t.star)[1]
#                     s3 <- sqrt(colVars(t.star, na.rm=TRUE))^3
#                     a = (m3/s3)/6
#                     w = qnorm(colSums(ifelse(t.star <= t.hat.mat, 1, 0),na.rm=TRUE)/(R+1))
#                     z.tilde.alpha.low = w + qnorm(alpha)
#                     alpha.tilde.low = pnorm(w + z.tilde.alpha.low/(1 - a*z.tilde.alpha.low))
#                     z.tilde.alpha.high= w + qnorm(1-alpha)
#                     alpha.tilde.high = pnorm(w + z.tilde.alpha.high/(1 - a*z.tilde.alpha.high))
#                     BCaL <- rbind(BCaL, interpolation(R=R, alpha=alpha.tilde.low, sort.mat=t.star.sort))
#                     BCaU <- rbind(BCaU, interpolation(R=R, alpha=alpha.tilde.high, sort.mat=t.star.sort))
#                     ## CI basic
#                     basicU <- rbind(basicU, 2*t.hat.s - interpolation(R=R, alpha=alpha, sort.mat=t.star.sort))
#                     basicL <- rbind(basicL, 2*t.hat.s - interpolation(R=R, alpha=1-alpha, sort.mat=t.star.sort))
#                     ##_________
# } #_____________________________ cond R
#                 if(!is.null(t.rea.s)) t.rea <- rbind(t.rea, t.rea.s)
#                 t.hat <- rbind(t.hat, t.hat.s)
#                 N <- list(N, list(N.s))
 
#                 sgreen <- rgb(0.067,0.431,0.22);    bgreen <- rgb(0.85, 0.85, 0.85);    bbgreen <- rgb(0.93,0.93, 0.93);    rolex <- rgb(0,0.25,0.067);    rohgreen <- rgb(0.46,0.80,0.30);    rolgreen <- rgb(0.66,1,0.46);    rohblue <- rgb(0,0.33,0.33);    roblue <- rgb(0.33,0.66,0.66);    vio <-rgb(167/256,0,0.5);    dgelb <- rgb(0.9,0.9,0.75);    nblue <- rgb(0,0,0.4);    nred <- rgb(0.4,0,0);    nred <- rgb(1,0,0.145);    bred <- rgb(1,0.5,0.5);    red <- rgb(1,0,0);    blue <- rgb(0,0,1);    grey <- rgb(0.5,0.5,0.5);    bgrey <- rgb(0.8,0.8,0.8);    ngrey <- rgb(0.2,0.2,0.2);    yellow <- rgb(0.9, 0.8, 0.5);    green <- rgb(0.2,0.5,0.2);    bg_red <- rgb(1,0.8,0.8);   bg_grey <- rgb(0.8,0.8,0.8)
# #################################################################################################### dPLOT
#                 if(any(output=="dPLOT") && s<5 && R>1){ # dPLOT cond
#                   la <- dim(percL)[1] # to select the last row
# ### BIAS.pdf
#                   pdf(paste("~/Desktop/",stat,"S",s,"_BIAS.pdf", sep=""), width=8, height=6)
#                   for (k in 1:dim(t.hat)[2]){
#                     par(mfrow=c(1,1), pty='s')
#                     plot(bias, ylim=c(0,22), xlim=c(0,60), xlab="", ylab="", main="", type="n", axes=FALSE)
#                     if (stat=="HER" | stat=="HER_ind") text(y=22, x=1, label=paste(dimnames(theta)[[1]]), pos=4,col=red)
#                     text(y=21, x=1, label=paste("trait:",(dimnames(theta)[[2]])[k]), pos=4,col=red)
#                     text(y=19, x=1, label=paste("R=",R, "Q=",mean(N.s,na.rm=TRUE)), pos=4,col=green)
#                     text(y=18, x=1, label=paste("transformation:",transformation), pos=4,col=green)
#                     text(y=16, x=1, label=paste("t:         ", round(utnf(theta[k]),3)), pos=4,col=green)
#                     if(!is.null(t.rea.s))text(y=15, x=1, label=paste("t.realized:", round(utnf(t.rea.s[k]),3)), pos=4,col=green)
#                     text(y=14, x=1, label=paste("t.hat:     ", round(utnf(t.hat.s[k]),3)), pos=4,col=green)
#                     text(y=11, x=1, label=paste(round(1-alpha,2)," - Confidence Intervals:"), pos=4,col=green)
#                     ## Sen
#                     if(!is.null(senL[la,k])&!is.null(senU[la,k])) text(y=9, x=1, label=paste("Sen92:  (", round(utnf(senL[la,k]),3), " ; ", round(utnf(senU[la,k]),3),")"), pos=4,col=green)
#                     ## perc
#                     if (!is.na(percL[la,k])&!is.na(percU[la,k])) text(y=8, x=1, label=paste("percentile:  (", round(utnf(percL[la,k]),3), " ; ", round(utnf(percU[la,k]),3),")"), pos=4,col=green)
#                     ## BCa
#                     text(y=7, x=1, label=paste("BCa:  \\ci{",if (!is.na(BCaL[la,k])){round(utnf(BCaL[la,k]),2)}else{NA},"}{",if(!is.na(BCaU[la,k])){round(utnf(BCaU[la,k]),2)}else{NA},"}"), pos=4,col=green)
#                     ## basic
#                     if (!is.na(basicL[la,k])&!is.na(basicU[la,k])) text(y=6, x=1, label=paste("basic:  \\ci{", round(utnf(basicL[la,k]),2), "}{", round(utnf(basicU[la,k]),2),"}"), pos=4,col=green)
#                     ## stud
#                     if(!is.na(studL[la,k])&!is.na(studU[la,k])) text(y=5, x=1, label=paste("studentized:  (", round(utnf(studL[la,k]),3), " ; ", round(utnf(studU[la,k]),3),")"), pos=4,col=green)
#                     #text(y=1, x=1, label=paste("numbers: sires:",design[4], "dams:",design[5], "ind:", design[6]), pos=4,col=green)
# ###
#                     truehist(bias[,k], xlab="bias: t.star - t.hat", main="bias")
#                     lines(density(bias[,k], na.rm=TRUE))
#                     abline(v=0)
# ###
#                     qqnorm(bias[,k], main="bias (normalQQ)")
#                     qqline(bias[,k])
# ###
#                     if(length(v.star[,k])>5){
#                       plot(t.star[,k], v.star[,k], type="n", main="Variance stable?")
#                       points(t.star[,k], v.star[,k], pch=paste(N.s))}
#                   }
#                   dev.off()
#                 } # dPLOT cond
#     } # S cond
#   } #S loop
# ###___________________________  S loop
# ############################################################################################# sPLOT
#   if(S>3) {
#     t.hat[t.hat<0] <- 0 # because otherwise we can not compare the mean and variance beteen REML and ANOVA
#     t.hat[t.hat>0.99] <- 0.99 # because heritability can not be larger than 1
#     bias.hat <- colMeans(t.hat - matrix(theta, nrow=dim(t.hat)[1],ncol=dim(t.hat)[2], byrow=TRUE)) #the bootstrap bias of the estimator
# #      bias.hat <- colMeans(t.hat - matrix(colMeans(t.hat), nrow=dim(t.hat)[1],ncol=dim(t.hat)[2], byrow=TRUE)) #the bootstrap bias of the estimator
#       theta <- matrix(colMeans(t.hat), nrow=1, ncol=1) # damit der BIAS nicht die error rates beeintraechtigt!
#       var.hat <- colVars(t.hat, na.rm=TRUE) #the bootstrap variance of the estimator
#       var.rea <- colVars(t.rea, na.rm=TRUE) #the bootstrap variance of the realization
#       mse <- bias.hat^2+var.hat 
# ### CI quality
#       perc <- CIquality(percL, percU, theta)
#       BCa <- CIquality(BCaL, BCaU, theta)
#       basic <- CIquality(basicL, basicU, theta)
#       stud <- CIquality(studL, studU, theta)
#       sen <- CIquality(senL, senU, theta)
#       for (k in 1:length(theta)) { # therefore the number of traits
#         if(sum(!is.na(t.hat[1:dim(t.hat)[1],k]))>1){den.t.hat <- density(t.hat[,k], na.rm=TRUE)}else{den.t.hat <- NULL}
#         if(sum(!is.na(t.rea[1:dim(t.rea)[1],k]))>1){den.t.rea <- density(t.rea[,k], na.rm=TRUE)}else{den.t.rea <- NULL}
#                                         #      xlimes <- c(min(c(den.t.rea$x,den.t.hat$x)),max(c(den.t.rea$x,den.t.hat$x)))
#                                         #      ylimes <- c(0,max(c(den.t.rea$y,den.t.hat$y)))
#                                         #ylimes <- c(0,max(c(den.t.rea$y,den.t.hat$y))) if est and rea should be shown side by side
#             ylimes <- NULL
#             #ylimes <- c(0,3)
# ### TEX file ###
#         ## change the name of the folder after changing the argument in the statistic()!
#         if(any(output=="tex")){
#         filenameCAT <- paste(path,"III_results","_",frommethod,"_",100*alpha,".tex",sep="")
#         cat(paste(strsplit(folder,split="--")[[1]][1],"&",strsplit(folder,split="--")[[1]][2],"&",format(round(utnf(theta[k]),3), nsmall=3),"&",
#         format(round(utnf(bias.hat[k]),3), nsmall=3),"&",
#         format(round(utnf(var.hat[k]),3), nsmall=3),"&",
#         format(round(utnf(var.rea[k]),3), nsmall=3),"&",
#         format(round(mean(utnf(senL[,k])),3), nsmall=3), "&",
#         format(round(mean(utnf(senU[,k])),3), nsmall=3), "&",
#         format(round(mean(utnf(senU[,k])-utnf(senL[,k])),3),
#         nsmall=3), "&", format(round(utnf(sen$Lerror[k]),3),
#         nsmall=3),"&", format(round(utnf(sen$Uerror[k]),3),
#         nsmall=3),"&", format(round(utnf(sen$error[k]),3),
#         nsmall=3),"\\NN%",folder,"\n"), file=filenameCAT, append=TRUE)
#       }#___tex
# ### CIlength
#         ##
#         if(any(output=="CIlength")){
# 	    CIlength <- list(utnf=utnf, which=folder,  frommethod=frommethod, alpha=alpha, boxDATA=boxplot(utnf(senU[,k]) - utnf(senL[,k])),  minLength=min(utnf(senU[,k])-utnf(senL[,k])),  maxLength=max(utnf(senU[,k])-utnf(senL[,k])), summaryLength=summary(utnf(senU[,k])-utnf(senL[,k])))
#         }
# ### sPLOTcolor ###
#         if(any(output=="sPLOTdata")){
#             pdf(paste("~/Desktop/","T",k,"_DIST",alpha,".pdf", sep=""), width=8, height=10)
#             par(mfrow=c(1,1), pty='m',bg=rgb(1,1,0.8),fg=ngrey,col.axis=ngrey, col.main=ngrey, mar=c(1,1,1,1))
#             ##
#             plot(t.hat, ylim=c(0,28), xlim=c(0,60), xlab="", ylab="", main="", type="n", axes=FALSE)
#             text(y=28, x=1, label=paste(dimnames(t.hat)[[1]])[k], pos=4,col=red)
#             text(y=27, x=1, label=paste("trait:",(dimnames(t.hat)[[2]])[k]), pos=4,col=red)
#             text(y=26, x=1, label=paste("&",round(utnf(theta[k]),3),"&", round(utnf(bias.hat[k]),3),"&", round(utnf(var.hat[k]),3),"&", round(utnf(var.rea[k]),3),"&", round(mean(utnf(senL[,k])),3), "&", round(mean(utnf(senU[,k])),3), "&", round(mean(utnf(senU[,k])-utnf(senL[,k])),3), "&", round(utnf(sen$Lerror[k]),3),"&", round(utnf(sen$Uerror[k]),3),"&", round(utnf(sen$error[k]),3)), pos=4, col=red)
#             text(y=25, x=1, label=paste("design:"), pos=4,col=ngrey)
#             text(y=24, x=1, label=paste("missing:",missing), pos=4,col=ngrey)
# #            text(y=23, x=1, label=paste("sires:",design[4], ",  dams:",design[5], ",  ind:", design[6]), pos=4,col=ngrey) 
#             text(y=21, x=1, label=paste("S=",S,"R=",R, "Q=",Q), pos=4,col=ngrey)
#             text(y=19, x=1, label=paste("transformation:",transformation), pos=4,col=ngrey)
#             ##
#             text(y=17, x=1, label=paste("bias / var / mse  ", round(utnf(bias.hat[k]),3)," / ", round(utnf(var.hat[k]),3)," / ", round(utnf(mse[k]),3)), pos=4,col=ngrey)
#             ##
#             text(y=16, x= 1, label=paste("Nominal Error Rate:", round(alpha,2)), pos=4,col=ngrey)
#             text(y= 15, x= 1, label=paste("Empirical Error Rates (mean; var)"), pos=4,col=ngrey)
#             text(y= 15, x=47, label=paste(" -> total"), pos=4,col=ngrey)
#             ## perc
#             if(!is.null(perc)){
#               text(y=14, x=1, label="percentile:", col=bred)
#               text(y=14, x=10, label=paste("lower:",round(utnf(perc$Lerror[k]),3), " ( ", round(mean(utnf(percL[,k]),na.rm=TRUE),3),";",round(var(utnf(percL[,k]),na.rm=TRUE),3),")"), pos=4, col=bred)
#               text(y=13, x=10, label=paste("upper:",round(utnf(perc$error[k]),3), " ( ", round(mean(utnf(percU[,k]),na.rm=TRUE),3),";",round(var(utnf(percU[,k]),na.rm=TRUE),3),")"), pos=4, col=bred)
#               text(y=13, x=47, label=paste(" -> ", round(utnf(perc$error[k]),3),""), pos=4,col=bred)
#             }
#             ## BCa
#             if(!is.null(BCa)){
#               text(y=12, x=1, label="BCa:", col=green)
#               text(y=12, x=10, label=paste("lower:",round(utnf(BCa$Lerror[k]),3), " ( ", round(mean(utnf(BCaL[,k]),na.rm=TRUE),3),";",round(var(utnf(BCaL[,k]),na.rm=TRUE),3),")"), pos=4, col=green)
#               text(y=11, x=10, label=paste("upper:",round(utnf(BCa$Uerror[k]),3), " \\ci{", round(mean(utnf(BCaU[,k]),na.rm=TRUE),3),"}{",round(var(utnf(BCaU[,k]),na.rm=TRUE),3),"}"), pos=4, col=green)
#               text(y=11, x=47, label=paste(" -> ", round(utnf(BCa$error[k]),3),""), pos=4,col=green)
#             }
#             ## basic
#             if(!is.null(basic)){
#               text(y=10, x=1, label="basic:", col=vio)
#               text(y=10, x=10, label=paste("lower:",round(utnf(basic$Lerror[k]),3), " ( ", round(mean(utnf(basicL[,k]),na.rm=TRUE),3),";",round(var(utnf(basicL[,k]),na.rm=TRUE),3),")"), pos=4, col=vio)
#               text(y=9, x=10, label=paste("upper:",round(utnf(basic$Uerror[k]),3), " ( ", round(mean(utnf(basicU[,k]),na.rm=TRUE),3),";",round(var(utnf(basicU[,k]),na.rm=TRUE),3),")"), pos=4, col=vio)
#               text(y=9, x=47, label=paste(" -> ", round(utnf(basic$error[k]),3),""), pos=4,col=vio)
#             }
#             ## stud
#             if(!is.null(stud)){
#               text(y=8, x=2, label="studentized:", col=red)
#               text(y=8, x=10, label=paste("lower:",round(utnf(stud$Lerror[k]),3), " ( ", round(mean(utnf(studL[,k]),na.rm=TRUE),3),";",round(var(utnf(studL[,k]),na.rm=TRUE),3),")"), pos=4, col=red)
#               text(y=7, x=10, label=paste("upper:",round(utnf(stud$Uerror[k]),3), " ( ", round(mean(utnf(studU[,k]),na.rm=TRUE),3),";",round(var(utnf(studU[,k]),na.rm=TRUE),3),")"), pos=4, col=red)
#               text(y=7, x=47, label=paste(" -> ", round(utnf(stud$error[k]),3),""), pos=4,col=red)
#             }
#             ## sen
#             if(!is.null(sen)){
#               text(y=6, x=1, label="Sen92:", col=blue)
#               text(y=6, x=10, label=paste("lower:",round(utnf(sen$Lerror[k]),3), " ( ", round(mean(utnf(senL[,k]),na.rm=TRUE),3),";",round(var(utnf(senL[,k]),na.rm=TRUE),3),")"), pos=4, col=blue)
#               text(y=5, x=10, label=paste("upper:",round(utnf(sen$Uerror[k]),3), " ( ", round(mean(utnf(senU[,k]),na.rm=TRUE),3),";",round(var(utnf(senU[,k]),na.rm=TRUE),3),")"), pos=4, col=blue)
#               text(y=5, x=47, label=paste(" -> ", round(utnf(sen$error[k]),3),""), pos=4,col=blue)
#             }
#             dev.off()
#           }#___sPLOTcolor
# ### DISTrea.pdf
#           #   if(all(missing==0)){ #with missing values t.rea has no meaning
# #               pdf(paste("~/Desktop/",stat,"_T",k,"_DISTrea.pdf", sep=""), width=4, height=3)
# #               par(mfrow=c(1,1), pty='s',bg=rgb(1,1,0.8),fg=ngrey,col.axis=ngrey, col.main=ngrey, mar=c(3,2,2,1))
# #               truehist(t.rea[,k], xlab="", ylab="", main="REALIZATIONS", col=bg_grey, ylim=ylimes, xlim=xlimes, axes=FALSE)
# #               mtext(statname, side=1, line=2, cex=1)
# #               mtext("Density", side=2, line=2, cex=1)
# #               if(class(den.t.rea)=="density") lines(den.t.rea, col=grey)
# #               abline(v=theta[k])
# #               axis(side=1, at=c(0,0.2,0.4,0.6,0.8,1))
# #               axis(side=1, at=c(seq(from=0, to=1, by=0.2)), labels=FALSE)
# #               axis(side=2)
# #               dev.off()
# #             }
# ### sPLOTbw ###
#         if(any(output=="sPLOTbw")){
#           bred <- rgb(1,1,1); green <- rgb(1,1,1); vio <- rgb(1,1,1); red <- rgb(1,1,1); blue <- rgb(0,0,0)
#           pdf(paste("~/Desktop/","T",k,"_DISTest.pdf", sep=""), width=9, height=6) # print
#           par(mfrow=c(1,2), pty='s',bg=rgb(1,1,1),fg=rgb(0,0,0),col.axis=rgb(0,0,0), col.main=rgb(0,0,0), mar=c(3,2,2,1)) # print
#           truehist(t.hat[,k], xlab="", ylab="", main="",col=rgb(0.7,0.7,0.7), ylim=ylimes, xlim=limes, axes=FALSE)
#           #mtext(st side=1, line=2, cex=1) 
#           mtext("Density", side=2, line=2, cex=1)
#           if(class(den.t.hat)=="density") lines(den.t.hat ,col=rgb(0,0,0))
#           if(est.ci){
#             if(class(perc$Udensity[[k]])=="density") lines(perc$Udensity[[k]], col=bred, lty=2)
#             if(class(perc$Ldensity[[k]])=="density") lines(perc$Ldensity[[k]], col=bred, lty=2)
#             if(class(BCa$Ldensity[[k]])=="density") lines(BCa$Ldensity[[k]], col=green, lty=2)
#             if(class(BCa$Udensity[[k]])=="density") lines(BCa$Udensity[[k]], col=green, lty=2)
#             if(class(basic$Ldensity[[k]])=="density") lines(basic$Ldensity[[k]], col=vio, lty=2)
#             if(class(basic$Udensity[[k]])=="density") lines(basic$Udensity[[k]], col=vio, lty=2)
#             if(class(stud$Ldensity[[k]])=="density") lines(stud$Ldensity[[k]], col=red, lty=2)
#             if(class(stud$Udensity[[k]])=="density") lines(stud$Udensity[[k]], col=red, lty=2)
#             if(class(sen$Ldensity[[k]])=="density") lines(sen$Ldensity[[k]], col=blue, lty=2)
#             if(class(sen$Udensity[[k]])=="density") lines(sen$Udensity[[k]], col=blue, lty=3, lwd=2.5)
#           }
#           abline(v=theta[k])
#             #text(x=xlimes[2], y=1*ylimes[2],label=paste("missing:"), pos=2)
#                                         #text(x=xlimes[2], y=0.9*ylimes[2],label=paste(round(missing,2)),pos=2)
#           axis(side=2)
#           axis(side=1, at=c(0,0.2,0.4,0.6,0.8,1.0))
#           axis(side=1, at=c(seq(from=0, to=1, by=0.2)), labels=FALSE)
#           qqnorm(t.hat[,k])
#           dev.off()
#         }
# ### sPLOTcolor
#         if(any(output=="sPLOTcolor")){
#           pdf(paste("~/Desktop/","T",k,"_DISTest",alpha,".pdf", sep=""), width=4, height=3) # presentation
#           par(mfrow=c(1,1), pty='s',bg=rgb(1,1,0.8),fg=ngrey,col.axis=ngrey, col.main=ngrey, mar=c(3,2,2,1)) #presentation
#           truehist(t.hat[,k], xlab="", ylab="", main="Estimation",col=bg_red, ylim=ylimes, xlim=statlimes, axes=FALSE)
#           #mtext(statname, side=1, line=2, cex=1) 
#           mtext("Density", side=2, line=2, cex=1)
#           if(class(den.t.hat)=="density") lines(den.t.hat ,col=bg_red)
#           if(est.ci){
#             if(class(perc$Udensity[[k]])=="density") lines(perc$Udensity[[k]], col=bred, lty=2)
#             if(class(perc$Ldensity[[k]])=="density") lines(perc$Ldensity[[k]], col=bred, lty=2)
#             if(class(BCa$Ldensity[[k]])=="density") lines(BCa$Ldensity[[k]], col=green, lty=2)
#             if(class(BCa$Udensity[[k]])=="density") lines(BCa$Udensity[[k]], col=green, lty=2)
#             if(class(basic$Ldensity[[k]])=="density") lines(basic$Ldensity[[k]], col=vio, lty=2)
#             if(class(basic$Udensity[[k]])=="density") lines(basic$Udensity[[k]], col=vio, lty=2)
#             if(class(stud$Ldensity[[k]])=="density") lines(stud$Ldensity[[k]], col=red, lty=2)
#             if(class(stud$Udensity[[k]])=="density") lines(stud$Udensity[[k]], col=red, lty=2)
#             if(class(sen$Ldensity[[k]])=="density") lines(sen$Ldensity[[k]], col=blue, lty=2)
#             if(class(sen$Udensity[[k]])=="density") lines(sen$Udensity[[k]], col=blue, lty=2)
#           }
#           abline(v=theta[k])
#                                         #text(x=xlimes[2], y=1*ylimes[2],label=paste("missing:"), pos=2)
#                                         #text(x=xlimes[2], y=0.9*ylimes[2],label=paste(round(missing,2)),pos=2)
#           axis(side=2)
#           axis(side=1, at=c(0,0.2,0.4,0.6,0.8,1.0))
#           axis(side=1, at=c(seq(from=0, to=1, by=0.2)), labels=FALSE)
#           dev.off()
#         }#___sPLOTbw
#         ###sPLOTphd
#         if(any(output=="sPLOTphd")){
#           pdf(file="/Users/thomasfabbro/Documents/Phd/1.Presentations/06.05.11_PhDtest/fig/DISest.pdf", width=pdfwidth2, height=pdfheight2)
#           par(bg=bgRGB, fg=fgRGB, col.axis=fgRGB, col.lab=fgRGB, col.main=fgRGB, mgp=c(2,1,0), mar=c(3,3,0.5,1))
#           hist(t.hat[,k], freq=FALSE, axes=FALSE,
#                xlab="Heritability", ylab="Density",
#                xlim=c(0,1), ylim=c(0,4.5),
#                col=col.1, main="")
#           axis(side=1)
#           axis(side=2)
#           abline(v=theta[k])
#           dev.off()
#           pdf(file="/Users/thomasfabbro/Documents/Phd/1.Presentations/06.05.11_PhDtest/fig/DISrea.pdf", width=pdfwidth2, height=pdfheight2)
#           par(bg=bgRGB, fg=fgRGB, col.axis=fgRGB, col.lab=fgRGB, col.main=fgRGB, mgp=c(2,1,0), mar=c(3,3,0.5,1))
#           hist(t.rea[,k], freq=FALSE, axes=FALSE,
#                xlab="Heritability", ylab="Density",
#                xlim=c(0,1), ylim=c(0,4.5),
#                col=col.100, main="")
#           axis(side=1)
#           axis(side=2)
#           abline(v=theta[k])
#           dev.off()
#         }
#       }
#     }
# ###OUTPUT
# if(any(output=="CIlength")){CIlength}
# }
