pcorr <- function(corr){
    if(all(!is.na(corr))){
	    corr <- - matrixscale(solve(corr))
	} else {
	    corr <- NA
	}
	corr	
}

# This function calculates the partial correlation matrix from a correlation matrix


## Test with
## Math Marks Data (available e.g. from Package ggm)
## > C <- cor(marks)
## > PCorr(C)
##                mechanics     vectors    algebra      analysis  statistics
## mechanics  -1.0000000000  0.32846282  0.2292442 -0.0007121818  0.02550897
## vectors     0.3284628169 -1.00000000  0.2816015  0.0778303605  0.01996942
## algebra     0.2292441885  0.28160149 -1.0000000  0.4317713938  0.35673607
## analysis   -0.0007121818  0.07783036  0.4317714 -1.0000000000  0.25277656
## statistics  0.0255089719  0.01996942  0.3567361  0.2527765574 -1.00000000


## result with the parcor-function from package ggm:
## > S <- var(marks)
## > parcor(S)
##                mechanics    vectors   algebra      analysis statistics
## mechanics   1.0000000000 0.32846282 0.2292442 -0.0007121818 0.02550897
## vectors     0.3284628169 1.00000000 0.2816015  0.0778303605 0.01996942
## algebra     0.2292441885 0.28160149 1.0000000  0.4317713938 0.35673607
## analysis   -0.0007121818 0.07783036 0.4317714  1.0000000000 0.25277656
## statistics  0.0255089719 0.01996942 0.3567361  0.2527765574 1.00000000



