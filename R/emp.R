emp <-function(data.use=NULL, fixedblock.use="fixedblock", character.use, environment.use="environment", randomblock.use="randomblock", sire.use="sire", dam.use="dam", individual.use="individual", without="", partitioning="REML", file=TRUE, path="~/Desktop/myproject.qgen/"){
### the arguments:         data.use=NULL; environment.use="environment";  sire.use="sire"; dam.use="dam"; individual.use="individual"; without=""; partitioning="REML"; file=TRUE; path="~/Desktop/myproject.Qgen/"
###
### default dataset: SINAPIS
### without a selected dataframe, the data from the included dataframe (called "sinapis") is used
  if(is.null(data.use)){
    data(sinapis)
    data.use <- sinapis
    fixedblock.use="fixedblock"
    character.use = c("leafarea","plantweight")
    envionment.use="environment"
    randomblock.use="randomblock"
    sire.use="sire"
    dam.use="dam"
    individual.use="individual"
  }
### DATA
  ## a dataframe is constructed where values from all traits are in one column
  mat <- data.use[,1] #adding the first column
  for (i in character.use){
    new <- data.frame(data.use[,fixedblock.use], rep(i,dim(data.use)[1]), data.use[,environment.use], data.use[,randomblock.use], data.use[,sire.use], data.use[,dam.use], data.use[,individual.use], data.use[,i])
    mat <- data.frame(mat, new)
  }
  matdat <- mat[,-1] #removing  the first column
  ## droping unused factor levels and rearranging
  data.out <- data.frame(fb=matdat[,1][drop=TRUE], ch=matdat[,2][drop=TRUE], en=matdat[,3][drop=TRUE], rb=matdat[,4][drop=TRUE], si=matdat[,5][drop=TRUE], da=matdat[,6][drop=TRUE], id=matdat[,7][drop=TRUE], y=matdat[,8][drop=TRUE])
  #removing all NA in y-coloumn
  data.out <- data.out[!is.na(data.out[,8]),] 
### DIAGNOSTIC PLOT
  ## if (length(character.use)>1){
  ## chix <- na.omit(data.frame(data.use[c(number,character.use)]))
  ##   chix <- chix[!(chix$number%in%without) ,]
  ## pdf("~/Desktop/chiplot.pdf")
  ##   Chiplot(chix, id=TRUE)
  ## dev.off()
  ## }
  ##
### SUPL information extraction
  chN <- length(levels(data.out$ch))
  enN <- length(levels(data.out$en))
  fbN <- length(levels(data.out$fb))
  siN <- length(levels(data.out$si)) #total number of sires
  daN <- max(tapply(data.out$da, data.out$si, function(x) length(tabulate(x)[tabulate(x)>0])),na.rm=TRUE)
  #
  data.out$si_da=as.factor(paste(data.out$si, data.out$da, sep="_"))
  idN <- max(tapply(data.out$id, list(data.out$si_da,data.out$en,data.out$fb), function(x) length(tabulate(x)[tabulate(x)>0])),na.rm=TRUE)
  data.out$si_da <- NULL
  ## calculates an array with the proportion of missing values
  miss <- 1-(tapply(data.out$y, list(data.out$en,data.out$ch, data.out$fb), length)/(siN*daN*idN))
### OUTPUT
  origEMP <- new("orig", hist="emp", warn="", time=date(), part=partitioning)
  suplEMP <- new("supl", chN=as.integer(chN), enN=as.integer(enN), fbN=as.integer(fbN), rbN=as.integer(1), siN=as.integer(siN), daN=as.integer(daN), idN=as.integer(idN), miss=miss)
  paraEMP <- new("para")
  DATAEMP <- new("DATA", dat=data.out[!(data.out$id%in%without),])
  specEMP <- new("spec")
  ##
  emp.paraDATA <- new("paraDATA", orig=origEMP, supl=suplEMP, para=paraEMP, DATA=DATAEMP, spec=specEMP)
  ## 
  if (file){
    dir.create(path=path, showWarnings = FALSE, recursive = TRUE)
    save(emp.paraDATA, file=paste(path, "emp.rda", sep=""))
  }
  emp.paraDATA
}
