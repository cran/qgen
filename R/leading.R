leading <- function(x){
  x <- as.integer(x)
  full <- paste("0000",x,sep="")#; print(full)
  out <- vector(mode="character", length=length(full))
  for (i in seq(along.with=full)){
    len <- length(strsplit(full,split="")[[i]])#; print(len)
    out[i] <- paste(substring(full[i], first=len-3, last=len))
  }
  return(out)
  }
## This function takes a numeric value or list of numeric values and creates a character string with a constant length of 4
## Examples:
## leading(4) -> "0004"
## leading(c(1,2,3,4))  returns  c("0001", "0002", "0003", "0004")

