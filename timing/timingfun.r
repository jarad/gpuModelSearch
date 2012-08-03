source("~/gpuModelSearch/Rwrapper/rlmsearch.r")
source("~/gpuModelSearch/Cwrapper/clmsearch.r")
source("~/gpuModelSearch/Csmart/cslmsearch.r")


##generate a model matrix
genX <- function(n, k){
  X <- matrix(c(rep(1,n),rnorm(n*k)),ncol=k+1)
  colnames(X) <- c("Int", paste("x", c(1:k), sep=""))
  return(X)
}

##generate a response vector
genY <- function(n,k,X){
  betas <- t(t(c(rnorm(k+1))))
  Y <- X%*%betas + c(rnorm(n, sd=5))
  return(Y)
}

gpuReset <- function(){

  if( !is.loaded("gpuReset") )
    dyn.load("devreset.so")

  out <- 0L

  id <- getGpuId()

  z <- .C("gpuReset", out=out)

  if(out == 0)
    out <- paste(c("GPU", id, "reset successfully"), collapse=" ")
  else
    out <- paste(c("GPU", id, "failed to reset"), collapse=" ")
  
  return(out)
}

