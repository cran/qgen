cal <- function(filename, repetitions=2, spliting=1, path="~/Desktop/myproject.qgen/") {
  dir.create(path=path, showWarnings = FALSE, recursive = TRUE)
### which file to create
  outputCAL <- paste(path, "CALfile.r", sep="")
### which group
    group <- "group1"
  if (substring(filename, first=7, last=7)=="S" & substring(filename, first=17, last=17)!="R"){group <- "group2"}  ## here all "multiTSyyyy.rda"
  if (substring(filename, first=17, last=17)=="R"){group <- "group3";print(filename)} ## here all multiTxxxxoyyyyRzzzz.rda"
  if (substring(filename, first=27, last=27)=="Q"){stop("? Do you really want to perform a tripple bootstrap ?\n-> this is not supportet; yet?!")}
### GROUP2
  if(group=="group2"){
    level <- "R"
    total.file.number <- 1
    total.run.number <- eval(parse(text=substring(filename, first=8, last=12)))
    cat("input.name <- load(\"",path,filename,"\")","\n", sep="",  file=outputCAL, append=FALSE)
    cat("eval(parse(text=paste(\"input\", \" <- \", input.name))); rm(list=c(input.name)); rm(input.name)","\n", sep="", file=outputCAL, append=TRUE)
    input.name <- load(file= paste(path, filename, sep=""))
    eval(parse(text=paste("input", " <- ", input.name)))
    ##
    for(run.number in 1:total.run.number){
      cat("add.multi <- vector(mode=\"list\", length=",repetitions,"); last.warning <- character()","\n", sep="",  file=outputCAL, append=TRUE)
      ##
      for(r in 1:repetitions){
        ##
        cat("add.paraDATA <- new(\"paraDATA\")", "\n", sep="",  file=outputCAL, append=TRUE)
        cat("add.paraDATA <- est(Sim(input@list.paraDATA[[",run.number,"]]))", "\n", sep="",  file=outputCAL, append=TRUE)
        cat("c(add.paraDATA@orig@warn, last.warning); c(add.paraDATA@orig@time, date()); ", "\n", sep="",  file=outputCAL, append=TRUE)
        cat("add.multi[[",r,"]] <- add.paraDATA", "\n", sep="",  file=outputCAL, append=TRUE)
        ##
      }
      ##
      cat("multi <- new(\"multi\", list.paraDATA=add.multi, level=\"",level,"\", x=",run.number,",y=",total.run.number,")","\n",  sep="", file=outputCAL,append=TRUE)
      cat(paste("save(multi, file=\"",path, "multiTS", leading(run.number), "o", leading(total.run.number), "R",leading(repetitions), ".rda\")", "\n", sep=""), sep="", file=outputCAL, append=TRUE)
      ##
    }
    cat("rm(add.paraDATA); rm(last.warning); rm(add.multi); rm(multi);rm(input)","\n",  sep="", file=outputCAL,append=TRUE)
  }
### GROUP3
  if(group=="group3"){
    #print("we are in group 3")
    level <- "Q"
    total.file.number <- eval(parse(text=substring(filename, first=13, last=16))); print(total.file.number)
    cat("\n", sep="",  file=outputCAL, append=FALSE) ##opens a new file
    for(file.number in 1:total.file.number){
      cat("input.name <- load(\"",paste(path, "multiTS", leading(file.number), "o", leading(total.file.number),"R", substring(filename, first=18, last=21), ".rda", sep=""),"\")", "\n", sep="",  file=outputCAL, append=TRUE)
                                        #substring(filename, first=1, last=regexpr("[0-9]{4}o[0-9]{4}R[0-9]{4}.rda", filename)-1), leading(file.number),"o",leading(total.file.number),".rda","\")","\n", sep="",  file=outputCAL, append=FALSE)
      cat("eval(parse(text=paste(\"input\", \" <- \", input.name))); rm(list=c(input.name)); rm(input.name)","\n", sep="", file=outputCAL, append=TRUE)
      input.name <- load(paste(path, "multiTS", leading(file.number), "o", leading(total.file.number),"R", substring(filename, first=18, last=21), ".rda", sep=""))
      eval(parse(text=paste("input", " <- ", input.name)))
      ## read from the object
      total.run.number <- length(multi@list.paraDATA)
      ##
      for(run.number in 1:total.run.number){
        cat("add.multi <- vector(mode=\"list\", length=",repetitions,"); last.warning <- character()","\n", sep="",  file=outputCAL, append=TRUE)
        ##
        for(r in 1:repetitions){
          ##
          cat("add.paraDATA <- new(\"paraDATA\")", "\n", sep="",  file=outputCAL, append=TRUE)
          cat("add.paraDATA <- est(Sim(input@list.paraDATA[[",run.number,"]]))", "\n", sep="",  file=outputCAL, append=TRUE)
          cat("c(add.paraDATA@orig@warn, last.warning); c(add.paraDATA@orig@time, date()); ", "\n", sep="",  file=outputCAL, append=TRUE)
          cat("add.multi[[",r,"]] <- add.paraDATA", "\n", sep="",  file=outputCAL, append=TRUE)
          ##
        }
        ##
        cat("multi <- new(\"multi\", list.paraDATA=add.multi, level=\"",level,"\", x=",run.number,",y=",total.run.number,")","\n",  sep="", file=outputCAL,append=TRUE)
        cat(paste("save(multi, file=\"",path, "multiTS", leading(file.number), "o", leading(total.file.number), "R", leading(run.number), "o", leading(total.run.number), "Q", leading(repetitions), ".rda\")", "\n", sep=""), sep="", file=outputCAL, append=TRUE)
        ##
      }
    }
    ##
    cat("rm(add.paraDATA); rm(last.warning); rm(add.multi); rm(multi);rm(input)","\n",  sep="", file=outputCAL,append=TRUE)
  }  
### GROUP 1
  if(group=="group1"){ ## here all without "multiTS..." name
    print(group)
    cat("input.name <- load(\"",path,filename,"\")","\n", sep="",  file=outputCAL, append=FALSE)
    cat("eval(parse(text=paste(\"input\", \" <- \", input.name))); rm(list=c(input.name)); rm(input.name)","\n", sep="", file=outputCAL, append=TRUE)
    ##
    input.name <- load(file= paste(path,filename, sep=""))
    eval(parse(text=paste("input", " <- ", input.name)))
    if(class(input)=="paraDATA") {
      hist <- input@orig@hist
      if(identical(hist, c("the", "sim", "est"))){stop("cal() was used with an argument with history: c(\"the\" \"sim\" \"est\"); use cal() on the object with history \"the\".")}
      if(identical(hist, "the")){
        level <-"S"
        ##
        #cat("input@orig@time <- date()","\n", sep="",  file=outputCAL, append=TRUE)
        cat("input.multi <- new(\"multi\", list.paraDATA=list(input), level=\"T\")","\n", sep="",  file=outputCAL, append=TRUE)
        cat("multi <- input.multi","\n", sep="",  file=outputCAL, append=TRUE)
        ##
        cat(paste("save(multi, file=\"",path,"multiT.rda\")","\n", sep=""),  file=outputCAL, append=TRUE)
        filename <- paste("multiTS",leading(repetitions),".rda", sep="")
      }
      if(identical(hist, c("emp", "est"))){
        level <- "R"
        ##
        #cat("input@orig@time <- date()","\n", sep="",  file=outputCAL, append=TRUE)
        cat("input.multi <- new(\"multi\", list.paraDATA=list(input), level=\"S\")","\n", sep="",  file=outputCAL, append=TRUE)
        cat("multi <- input.multi","\n", sep="",  file=outputCAL, append=TRUE)
        ##
        cat(paste("save(multi, file=\"",path,"multiTS0001.rda\")","\n", sep=""),  file=outputCAL, append=TRUE)
#        filename <- "multiTS0001.rda"
        filename <- paste("multiTS0001o0001R",leading(repetitions),".rda", sep="")
      }
    }
    cat("add.multi <- vector(mode=\"list\", length=",repetitions,"); last.warning <- character()","\n", sep="",  file=outputCAL, append=TRUE)
    for(r in 1:repetitions){
      ##
      cat("add.paraDATA <- new(\"paraDATA\")", "\n", sep="",  file=outputCAL, append=TRUE)
      cat("add.paraDATA <- est(sim(input.multi@list.paraDATA[[1]]))", "\n", sep="",  file=outputCAL, append=TRUE)
      cat("c(add.paraDATA@orig@warn, last.warning); c(add.paraDATA@orig@time, date()); ", "\n", sep="",  file=outputCAL, append=TRUE)
      cat("add.multi[[",r,"]] <- add.paraDATA", "\n", sep="",  file=outputCAL, append=TRUE)
      ##
    }
    ##
    cat("multi <- new(\"multi\", list.paraDATA=add.multi, level=\"",level,"\")","\n",  sep="", file=outputCAL,append=TRUE)
    cat("save(multi, file=\"",path, filename,"\")", "\n", sep="", file=outputCAL, append=TRUE)
    cat("rm(add.paraDATA); rm(last.warning); rm(add.multi); rm(multi);rm(input.multi)","\n",  sep="", file=outputCAL,append=TRUE)
  }
### INFO
  print(paste("now you can source() the file:", outputCAL))
}
