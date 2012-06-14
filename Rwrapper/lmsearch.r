
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





gpuLmsearch <- function(Y, X, g, sortby="MargLike", storemodels=FALSE){
  p <- ncol(X) #number of regressors, including intercept
  k <- p - 1 #number of possible covariates
  n <- nrow(X) #sample size
  M <- 2^k #number of possible models (all contain intercept)
  
  if(storemodels = TRUE) ##only used if full model info is desired
    models <- list()

  ##C is a common multiplicative constant to marginal likelihood of all models
  Ynorm <- sum( ( Y-mean(Y) )^2 )
  C <- margC(n, Ynorm)
  

  ##first, center X:
  if(p>1)
    for(i in 2:p){
      mn <- mean(X[,i])
      X[,i] <- X[,i] - mn
    }

  ##initialize variables
  ID <- 1:M 
  BinId <- paste(ID) 
  Aic <-  ID 
  Bic <- ID
  MargLike <- ID
  Vars <- paste(ID)

  ##for each model
  for(i in ID){
    binid <- modelid(i,k)  #find binary ID of model
    Xm <- modelmat(X, binid)  #find model matrix for model i
    km <- sum(binid)  #number of covariates in model i
    pm <- km+1  #num covariates + 1 for intercept, in model i

    gout <- gpuLm.fit(Xm, Y) #fit model i on gpu using gputools
    ssr <- sum(gout$residuals ^ 2) #sum of squared residuals
    sighat <- ssr/gout$df.residual #estimate of sigma^2
    Rsq <- 1 - ssr/Ynorm #R^2

    ##save full model info, if wanted
    if(storemodels = TRUE){
      gout$BinID <- binid
      gout$ID <- i 
      models[[i]] <- gout ##full model information for model i
    }

    ##save model selection info
    BinId[i] <- paste(binid, collapse="")
    Aic[i] <- aic(n,pm,sighat)
    Bic[i] <- bic(n,pm,sighat)
    MargLike[i] <- marglik(n, km, g, Rsq, C)
    Vars[i] <- paste(colnames(Xm),collapse=" ")
  }

  Ar <- rank(Aic) ##create rankings by AIC
  Br <- rank(Bic) ##create rankings by BIC
  Mr <- rank(-MargLike) ##create rankings by Marginal Likelihood

  ##create output list, out$list contains data frame of model selection info
  ##out$models contains model specific information - coefficients, resids, et
  out <- list()
  out$list <- data.frame(ID=ID, BinaryID=BinId, AIC=Aic, AICrank=Ar, BIC=Bic,
                         BICrank=Br, MargLike=MargLike, MLrank=Mr, Variables=Vars)
  if(storemodels = TRUE)
    out$models <- models

  if(sortby=="AIC")
    out$list <- out$list[ order(Ar),]
  else if(sortby=="BIC")
    out$list <- out$list[ order(Br),]
  else
    out$list <- out$list[ order(Mr),]
  
  return(out)
}

##note: AIC and BIC defined only up to additive constant
##calculates AIC, up to additive constant
aic <- function(n, p, sighat){
  out <- 2*(p+1) + n*log(sighat)
  return(out)
}

##calculates BIC, up to additive constant
bic <- function(n, p, sighat){
  out <- (p+1)*log(n)+n*log(sighat)
  return(out)
}

##Calculates the marginal likelihood for a given model
##Formula from Liang et al 2007, pg 6, assumes constant g
##later it might be worth implementing with priors chosen for g
marglik <- function(n, k, g, Rsq, C){
  ##if we have the null model, we know Rsq = 0
  if(k==0)
    Rsq <- 0
  
  out <- (1+g) ^ ( (n-1-k)/2 )
  out <- out / (   ( 1 + g*(1-Rsq) ) ^ ( (n-1)/2 ) )
  out <- C * out
  return(out)
}

##common multiplicative constant on all marginal likelihoods
margC <- function(n, Ynorm){
  Ynorm <- sqrt( Ynorm  )^(-n+1)
  K <- gamma((n-1)/2) / pi^((n-1)/2) / sqrt(n)
  out <- K*Ynorm
  return(out)
}
