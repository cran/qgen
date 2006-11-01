sta <- function(filename="", path="~/Desktop/myproject.qgen/", statistic.name,  alpha=0.05,  transformation="none"){
  ## depending on the length of the filename...
  file.name.length <- length(strsplit(filename, split="")[[1]]); print(paste("file name has length: ",file.name.length))
  if(file.name.length==10){
    total.file.number <- 1
      total.file.number.second <- 1
      rep.length <- 1
      file.for.parse <- "print(paste(path,\"multiT.rda\", sep=\"\"))"
      id.for.parse <- "print(paste(\"T1\", sep=\"\"))"
    }
    if(file.name.length==15){
      total.file.number <- 1
      total.file.number.second <- 1
      rep.length <-eval(parse(text=substring(filename, first=8, last=12)))
      file.for.parse <- "print(paste(path,\"multiTS\",leading(rep.length),\".rda\", sep=\"\"))"
      id.for.parse <- "print(paste(\"T1S\", leading(rep.number), sep=\"\"))"
    }
    if(file.name.length==25){
      total.file.number <- eval(parse(text=substring(filename, first=13, last=16)))
      total.file.number.second <- 1
      rep.length <-eval(parse(text=substring(filename, first=18, last=22)))
      file.for.parse <- "print(paste(path,\"multiTS\",leading(file.number), \"o\",leading(total.file.number),\"R\",leading(rep.length),\".rda\", sep=\"\"))"
      id.for.parse <- "print(paste(\"T1S\", leading(file.number),\"R\",leading(rep.number), sep=\"\"))"
    }
    if(file.name.length==35){
      total.file.number <- eval(parse(text=substring(filename, first=13, last=16)))
      total.file.number.second <-  eval(parse(text=substring(filename, first=23, last=26)))
      rep.length <-eval(parse(text=substring(filename, first=28, last=32)))
      file.for.parse <- "print(paste(path,\"multiTS\",leading(file.number), \"o\",leading(total.file.number),\"R\", leading(file.number.second), \"o\", leading(total.file.number.second), \"Q\", leading(rep.length),\".rda\", sep=\"\"))"
      id.for.parse <- "print(paste(\"T1S\", leading(file.number),\"R\", leading(file.number.second),\"Q\",leading(rep.number), sep=\"\"))"
    }
    ## transformations
    if (transformation=="none") {tnf <- function(x) x;   utnf <- tnf}
    if (transformation=="logodds") {tnf <- function(x) if (min(x)<0) rep(NA, length(x)) else  log((x)/(1-x));  utnf <- function(y) (2.718^y)/(1+2.718^y)}
    if (transformation=="asin") {tnf <- function(x) if (all(x<1)) asin(sqrt(x)) else rep(NA, length(x)); utnf <- function(y) (sin(y))^2}
    if (transformation=="fisher") {tnf <- function(x) if(all(x<1) & all(x>-1)) 0.5*log((1+x)/(1-x)); utnf <- function(y) (2.718^(2*y)-1)/(1+2.718^(2*y))}
    ##
    stat.matrix <- matrix(nrow=rep.length*total.file.number.second*total.file.number, ncol=attributes(statistic.name)$stat.dim)
    dim.stat.matrix <- vector(mode="character", length=rep.length*total.file.number*total.file.number.second)
    for (file.number in 1:total.file.number){
      for (file.number.second in 1:total.file.number.second){
        input.name <- load(file=eval(parse(text=file.for.parse)))
        eval(parse(text=paste("input", " <- ", input.name)))
        if(class(input)!="multi") stop("the file you specified does not contain an object of class \"multi\"")
        statname <- attributes(statistic.name)$name
        rm(list=c(input.name))
        rm(input.name)
        ##
        #print(paste("rep.length; ", rep.length))
        #print(paste("total.file.number; ", total.file.number))
        #print(paste("total.file.number.second; ", total.file.number.second))
        print(paste("dim of stat.matrix; ", dim(stat.matrix)))
        ##print(paste("file.number; ", file.number))
        ##
        for (rep.number in 1:rep.length){
          #print(paste("which slot; ",rep.number + (rep.length*(file.number.second-1)) + (rep.length*total.file.number.second*(file.number-1))        ))
             dim.stat.matrix[rep.number + (rep.length*(file.number.second-1)) + (rep.length*total.file.number.second*(file.number-1))] <- eval(parse(text=id.for.parse))
          stat.matrix[rep.number + (rep.length*(file.number.second-1)) + (rep.length*total.file.number.second*(file.number-1)),] <- tnf(statistic.name(input@list.paraDATA[[rep.number]], alpha=alpha)@stat)
          #print(tnf(statistic.name(input@list.paraDATA[[rep.number]], alpha=alpha, frommethod=frommethod)@stat))
           }
      }
    }
    dimnames(stat.matrix) <- list(dim.stat.matrix, statname)
  dir.create(path=path, showWarnings = FALSE, recursive = TRUE)
  save(stat.matrix,file=paste(path,"stat", input@level,".rda",sep=""))
}
