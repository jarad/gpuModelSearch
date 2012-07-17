library(gputools)

gpuclmsearch <- function(X, Y, g=NULL, sort="AIC", nsave=1000){

  if(is.null(g))
    g <- max((ncol(X)-1)^2, nrow(X))

  out <- gpuclmsearch.fit(X, Y, g, sort, nsave)
  return(out)
}

gpuclmsearch.fit <- function(X, Y, g, sorttype, nsave){

  if(!is.loaded("lmsearch"))
    dyn.load("lmsearch.so")
  nsave <- min(nsave, 2^(ncol(X)-1))
  mode(g) <- "integer"
  mode(X) <- "single"
  mode(Y) <- "single"
  mode(nsave) <- "integer"
  n <- as.integer(nrow(X))
  p <- as.integer(ncol(X))
  ycols <- as.integer(ncol(Y))

  if(sorttype=="AIC"){
    sort <- 1L
  }
  else if(sorttype=="BIC"){
    sort <- 2L
  }
  else{
    sort <- 3L
  }
  
  aics <- rep(1000000000000, nsave)
  bics <- aics
  lmls <- -aics
  models    <- integer(nsave)
  probs     <- single(nsave)
  otherprob <- single(1)
  mode(aics) <- "single"
  mode(bics) <- "single"
  mode(lmls) <- "single"

  z <- .C("lmsearch", X, n, p, Y, ycols, g, aic=aics, bic=bics, lml=lmls,
          probs=probs, otherprob=otherprob, id=models, nsave, sort)

  attr(z$aic,       "Csingle") <- NULL
  attr(z$bin,       "Csingle") <- NULL
  attr(z$lml,       "Csingle") <- NULL
  attr(z$probs,     "Csingle") <- NULL
  attr(z$otherprob, "Csingle") <- NULL

  binid <- modelid(z$id, p-1, nsave)
  
  
  out <- data.frame(ID=z$id,  BinaryID=modelidchar(binid, nsave),
                    AIC=z$aic,
                    AICrank=rank(z$aic,ties.method="first"),
                    BIC=z$bic,
                    BICrank=rank(z$bic,ties.method="first"),
                    LogMargLike=z$lml,
                    LMLrank=rank(-z$lml,ties.method="first"),
                    PostProb=z$probs,
                    Variables=modelvars(binid, colnames(X)[-1], nsave))
    
 
  if(sorttype=="AIC"){
    out <- out[order(out$AICrank),]
  }
  else if(sorttype=="BIC"){
    out <- out[order(out$BICrank),]
  }
  else{
    out <- out[order(out$LMLrank),]
  }

  out <- list(Models=out, OtherProb=z$otherprob)
  
  return(out)
}


modelidchar <- function(binid, nsave){
  out <- rep("a",nsave)
  for(i in 1:nsave){
    out[i] <- paste(binid[i,], collapse="")
  }
  return(out)
}


modelvars <- function(binid, colnam, nsave){
  out <- rep("a", nsave)
  
  for(i in 1:nsave){
    out[i] <- paste(colnam[which(binid[i,]==1)], collapse=" ")
  }
  return(out)
}


modelid <-  function(id, k, nsave){
  M <- 2^k
  binid <- matrix(0,ncol=k, nrow=nsave)

  for(i in 1:nsave){
    remain <- id[i] - 1
    for(j in 1:k){
      tmp <- 2^(k - j)
      if(remain >= tmp){
        remain <- remain - tmp
        binid[i,j] <- 1
      }
    }
  }
  return(binid)
}


##generate a model matrix
genX <- function(n, k){
  X <- matrix(c(rep(1,n),rnorm(n*k)),ncol=k+1)
  colnames(X) <- c("Int", paste("x", c(1:k), sep=""))
  return(X)
}


##generate a response vector
genY <- function(n,k,X){
  betas <- t(t(c(rnorm(k+1))))
  Y <- X%*%betas
  return(Y)
}

