library(gputools)

gpuClmsearch <- function(X, Y, g=NULL, sort="AIC", nsave=1000){

  if(is.null(g))
    g <- max((ncol(X)-1)^2, nrow(X))

  out <- gpuClmsearch.fit(X, Y, g, sort, nsave)
  return(out)
}

gpuClmsearch.fit <- function(X, Y, g, sorttype, nsave){

  if(!is.loaded("Clmsearch"))
    dyn.load("~/gpuModelSearch/Cwrapper/lmsearch.so")
  nsave <- min(nsave, 2^(ncol(X)-1))
  mode(g) <- "integer"
  mode(X) <- "single"
  mode(Y) <- "single"
  mode(nsave) <- "integer"
  n <- as.integer(nrow(X))
  p <- as.integer(ncol(X))
  ycols <- as.integer(ncol(Y))

  binids <- matrix(0L,ncol=p-1, nrow=nsave)

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
  probs     <- rep(0, nsave)
  otherprob <- 0
  mode(aics) <- "single"
  mode(bics) <- "single"
  mode(lmls) <- "single"

  z <- .C("Clmsearch", X, n, p, Y, ycols, g, aic=aics, bic=bics, lml=lmls,
          probs=probs, otherprob=otherprob, id=models, bin=binids, nsave, sort)

  attr(z$aic,       "Csingle") <- NULL
  attr(z$bin,       "Csingle") <- NULL
  attr(z$lml,       "Csingle") <- NULL
  attr(z$bin,       "Csingle") <- NULL


  out <- data.frame(ID=z$id,  BinaryID=modelidchar(z$bin, nsave),
                    AIC=z$aic,
                    AICrank=rank(z$aic,ties.method="first"),
                    BIC=z$bic,
                    BICrank=rank(z$bic,ties.method="first"),
                    LogMargLike=z$lml,
                    LMLrank=rank(-z$lml,ties.method="first"),
                    PostProb=z$probs,
                    Variables=modelvars(z$bin, colnames(X)[-1], nsave))
    
 
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
    out[i] <- paste(c("Int", colnam[which(rev(binid[i,])==1)]), collapse=" ")
  }
  return(out)
}


