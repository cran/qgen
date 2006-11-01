interpolation <- function(R, alpha, sort.mat){
  if(!is.matrix(sort.mat)){stop("sort.mat is not of class matrix")}
  if(alpha>1||alpha<0){stop("the argument \"alpha\" must be between zero and one ")}  
  interpol <- vector()
  alphaa <- matrix(alpha, nrow=1, ncol=dim(sort.mat)[2])
  for  (z  in 1:dim(sort.mat)[1]){
#    print(paste("z:",z)); print(alphaa[z])
    if(!is.nan(alphaa[z])){
      k <- trunc((R+1)*alphaa[z])#;print(paste("k:",k))
      sorti <- matrix(sort.mat[z,],nrow=dim(sort.mat)[2],ncol=1)#; print(sorti)
      if (k==0 || k==R || k==(R+1)) {
        interpol.new <- NA
      }else{
        interpol.new <- sorti[k,]+ ((qnorm(alphaa[z])-qnorm(k/(R+1)))/(qnorm((k+1)/(R+1))-qnorm(k/(R+1))))*(sorti[k+1,]-sorti[k,])
      }
    }else{
      interpol.new <- rep(NA, times=dim(sort.mat)[2])}
    interpol <- c(interpol, interpol.new)
  }
  interpol
}
# After Davison & Hinkley; Bootstrap methods and their application; p195; eq: 5.8

#sort.mat is an array where the first dimension is S and the second is R
#this function then sorts only within the S (which is the same as among the R's)

#alpha can be a single number or a vector
