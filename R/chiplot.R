chiplot <- function(y, id=FALSE ) {# this function generates a Chi-square probability plot
    if(id){x <- y[,-1]; name <- y[,1]}else{ x <- y; name <- NULL}
    n <- nrow(x)#number of replicates
    p <- ncol(x)#number of variables
    xbar <- apply(x, 2, mean)#vector of means
    S <- var(x)#covariance matrix
    index <- (1:n)/(n+1)#sample quantiles
    dist <- mahalanobis(x,xbar,S)#get Mahalanobis distances on the easy way
    quant <- qchisq(index,p)
    shift <- 0.1*max(quant)
    plot(quant, dist[order(dist)], ylab = "Ordered distances", xlab = "Chi-square quantile", lwd=2,  pch=3)
    text(quant-shift, dist[order(dist)], name[order(dist)],cex=0.7)
}
