#### R functions #####


invlogit<-function(x) {1/(1+exp(-(x)))} # inverse logit

logit<-function(x) {log(x/(1-x))}

CV <- function(x) (sd(x,na.rm=TRUE)/mean(x,na.rm=TRUE)) # Coeffvcient Variation

se <- function(x) sd(x,na.rm=TRUE) / sqrt(length(x)) # stabdard error
#The standard error is the standard deviation of the mean in repeated samples from a population.

#### make colors transparent
#note: always pass alpha on the 0-255 scale
makeTransparent<-function(someColor, alpha=100)
{
  newColor<-col2rgb(someColor)
  apply(newColor, 2, function(curcoldata){rgb(red=curcoldata[1], green=curcoldata[2],
                                              blue=curcoldata[3],alpha=alpha, maxColorValue=255)})
}
# example: makeTransparent(col[1],20)






## Convert arrays to list of matrix
conv <- function(x,j){
  y <- list()
  for (i in 1:sum(Ncond)){
    res <- matrix(NA,ncol=dim(x[[1]][[1]])[1],nrow=nSIMUL)
    for (repmod in 1:nSIMUL){
      res[repmod,] <- x[[i]][[repmod]][,j]
    }
    y[[i]]<-res
  }
  return(y)
}



# Measuring Response to Selection:
# * The selection differential S is the difference of the base population mean and the mean of the selected parents. 
# * The selection response R is how much gain you make when mating the selected parents.
# * R = h2 x S
# * Intensity of Selection: This is just a standardized version of the selection differential; i.e., the selection differential divided by the standard phenotypic deviation of the trait. -> i
# * how adding harvest to natural mortality affects Ne, census size (N), and the ratio Ne/N
# variables: c("Lf", "gPercF", "gSLmid", "galphaS", paste(c("gFmid", "pFmid"), rep(1:4, each = 2), sep = ""), "gNeutral")

# Metrics to measure evolution
Selection_differential <- function(x,y){ # Differential seelction
  return(mean(y)-mean(x))
}

R <- function(x,y,h2){ # Response to selection
  return((mean(y)-mean(x))*h2)
}

Selection_intensity <- function(x,y){ # Intensity of selection
  return((mean(y)-mean(x))/sd(x))
}


# haldane <- function(x,t1,t2){
#   return((log(x[t1])-log(x[t2]))/
#            }


