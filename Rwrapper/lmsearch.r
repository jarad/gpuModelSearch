
library(gputools)

##converts model id to a binary string
modelid <-  function(id, k){
  M <- 2^k
  binid <- rep(0,k)
  remain <- id - 1
  for(i in 1:k){
    tmp <- 2^(k - i)
    if(remain >= tmp){
      remain <- remain - tmp
      binid[i] <- 1
    }
  }
  return(binid)
}

##takes model id and full model matrix, outputs model matrix for model id
modelmat <- function(X, binid){
  idx <- c(1, (which(binid==1)+1) )
  Xm <- matrix(  X[,idx], ncol=sum(binid)+1)
  colnames(Xm) <- colnames(X)[idx]
  return(Xm)
}




##needs to be turned into a function
gpuLmsearch <- function(Y, X, g, sortby="MargLike"){
  p <- ncol(X)
  k <- p-1
  n <- nrow(X)
  M <- 2^k
  models <- list()
  Ynorm <- sum( ( Y-mean(Y) )^2 )
  C <- margC(n, Ynorm)

  ##first, center X:
  if(p>1)
    for(i in 2:p){
      mn <- mean(X[,i])
      X[,i] <- X[,i] - mn
    }
  
  ID <- 1:M
  BinId <- paste(ID)
  Aic <-  ID
  Bic <- ID
  MargLike <- ID
  Vars <- paste(ID)
  
  for(i in ID){
    binid <- modelid(i,k)
    Xm <- modelmat(X, binid)
    km <- sum(binid)
    pm <- km+1
    ##dat <- data.frame(Y,Xm)
    gout <- gpuLm.fit(Xm, Y)
    ssr <- sum(gout$residuals ^ 2)
    sighat <- ssr/gout$df.residual
    Rsq <- 1 - ssr/Ynorm
    gout$ID <- i
    gout$BinID <- binid
    models[[i]] <- gout
    BinId[i] <- paste(binid, collapse="")
    Aic[i] <- aic(n,pm,sighat)
    Bic[i] <- bic(n,pm,sighat)
    MargLike[i] <- marglik(n, km, g, Rsq, C)
    Vars[i] <- paste(colnames(Xm),collapse=" ")
  }

  Ar <- rank(Aic)
  Br <- rank(Bic)
  Mr <- rank(-MargLike)
  
  out <- list()
  out$list <- data.frame(ID=ID, BinaryID=BinId, AIC=Aic, AICrank=Ar, BIC=Bic,
                         BICrank=Br, MargLike=MargLike, MLrank=Mr, Variables=Vars)
  out$models <- models

  if(sortby=="AIC")
    out$list <- out$list[ order(Ar),]
  else if(sortby=="BIC")
    out$list <- out$list[ order(Br),]
  else
    out$list <- out$list[ order(Mr),]
  
  return(out)
}

#note: AIC and BIC defined only up to additive constant
aic <- function(n, p, sighat){
  out <- 2*(p+1) + n*log(sighat)
  return(out)
}

bic <- function(n, p, sighat){
  out <- (p+1)*log(n)+n*log(sighat)
  return(out)
}

##Formula from Liang et al 2007, pg 6, assumes constant g
##later it might be worth implementing with priors chosen for g
marglik <- function(n, k, g, Rsq, C){
  if(k==0)
    print(Rsq)
  
  out <- (1+g) ^ ( (n-1-k)/2 )
  out <- out / (   ( 1 + g*(1-Rsq) ) ^ ( (n-1)/2 ) )
  out <- C * out
  return(out)
}

margC <- function(n, Ynorm){
  Ynorm <- sqrt( Ynorm  )^(-n+1)
  K <- gamma((n-1)/2) / pi^((n-1)/2) / sqrt(n)
  out <- K*Ynorm
  return(out)
}
