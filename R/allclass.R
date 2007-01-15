.onLoad <- function(lib, pkg) { }

################################################## paraDATA
setClass("orig",
         representation(hist="character", warn="character", time="character", part="character"),
         validity = function(object) {
           if(length(object@hist)!=length(object@warn)) return("the vectors in slots hist and warn should have the same length")
           if(length(object@hist)!=length(object@time)) return("the vectors in slots hist and time should have the same length")
           if(length(object@warn)!=length(object@time)) return("the vectors in slots warn and time should have the same length")
           TRUE
         })
##
setClass("supl",
         representation(chN="integer", enN="integer", fbN="integer", rbN="integer", siN="integer", daN="integer", idN="integer", miss="array"),
         validity = function(object) {
           ## miss has three dimensions (fixedblock, environment, character)
           if (length(dim(object@miss))!=3){
             return("the slot \"miss\" of \"supl\" should have three dimentions: one for fixedblock, environment, and character")
           }
           ## the dimensions of miss must be chN, enN, fbN
##            if (dim(object@miss)[1]!=object@enN){ #& dim(object@miss)[2]== object@supl@enN & dim(object@miss)[3]==object@supl@fbN)){
##              return("the slot \"miss\" of \"supl\" should have the dimentions c(enN, chN, fbN)")
##            }
           ## there must be one character
           Dim <- object@chN
           if (Dim < 0){
             return("chN must be larger or equal to 1")
           }
           ## 'else'	ok :
           TRUE
	 })
       # setValidity("supl", valiDsupl)
##
setClass("DATA", representation(dat="data.frame"))
##
setClass("para",
         representation(rbS="matrix", siS="matrix", daS="matrix", idS="matrix", phS="matrix", error="numeric", fixe="array"),
         validity=function(object){
           if(!identical(dim(object@rbS), dim(object@siS)))
             return("missmach in dimention of Sigma")
           TRUE
         })
# the the slot fixe is also used to transport the names of the traits: "environments" and "characters"
##
setClass("spec", representation(additional.partitioning="list", unbalanced="list", modelsummary="list", secondcontrast="list"))
##
setClass("paraDATA",
         representation(orig="orig",supl="supl",para="para",DATA="DATA",spec="spec")
         )
################################################## multi
setClass("multi",
         representation(list.paraDATA="list", level="character", x="numeric", y="numeric"))
################################################## stat
setClass("stat",
         representation(orig="orig", stat="numeric", lower.ci="numeric", upper.ci="numeric", lower.limes="numeric", upper.limes="numeric"))








## new("paraDATA")
# #d1 <- new("paraDATA", hist="na", part="n1", chN=5, enN=6, rbS=matrix(c(4,4,4,4)), siS=matrix(c(3,3,3,3)), dat=data.frame(ch=runif(45)))
# #dd <- new("orig", hist="na", part="n1")
# #d2 <- new("paraDATA", oriG=dd, chN=5, enN=6, rbS=matrix(c(4,4,4,4)), siS=matrix(c(3,3,3,3)), dat=data.frame(ch=runif(45)))
# #
# #whatis <- function(object) paste("an object of class:", data.class(object))
# #
# #setMethod("whatis", "paraDATA",
# #         function(object){
# #           paste("an object of class:", data.class(object), "(paraDATA?)")
# #         }
# #       )
