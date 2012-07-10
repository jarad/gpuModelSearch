
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




##Stores model selection info of the best 1000 models according to sortby
gpurlmsearch <- function(X, Y, g=nrow(Y), sortby="AIC", storemodels=FALSE,
                        printi=FALSE){
  p <- ncol(X) #number of regressors, including intercept
  k <- p - 1 #number of possible covariates
  n <- nrow(X) #sample size
  M <- 2^k #number of possible models (all contain intercept)
  nlist <- min(M, 1000)
  
  if(storemodels == TRUE) ##only used if full model info is desired
    models <- list()

  Ynorm <- sum( ( Y-mean(Y) )^2 )
  ##C is a common multiplicative constant to marginal likelihood of all models
  ## C <- margC(n, Ynorm)
  

  ##first, center X:
  ##if(p>1)
   ## for(i in 2:p){
    ##  mn <- mean(X[,i])
     ## X[,i] <- X[,i] - mn
   ## }

  ##initialize variables
  i <- 1
  ID <- 1:nlist
  BinId <- paste(ID) 
  Aic <-  rep(Inf, nlist) 
  Bic <- rep(Inf, nlist)
  LogMargLike <- rep(-Inf, nlist)
  Vars <- paste(ID)
  maxlml <- -Inf


  ##for each model
  while(i <= M){
    binid <- modelid(i,k)  #find binary ID of model
    Xm <- modelmat(X, binid)  #find model matrix for model i
    km <- sum(binid)  #number of covariates in model i
    pm <- km+1  #num covariates + 1 for intercept, in model i

    gout <- gpuLm.fit(Xm, Y) #fit model i on gpu using gputools
    ssr <- sum(gout$residuals ^ 2) #sum of squared residuals
    sighat <- ssr/gout$df.residual #estimate of sigma^2
    Rsq <- 1 - ssr/Ynorm #R^2

    lml <- logmarglik(n, km, g, Rsq)
    a <- aic(n,pm,sighat)
    b <- bic(n,pm,sighat)

    ##save model selection info
    if(sortby=="AIC"){
      WorstIdx <- which.max(Aic)
      WorstScore <- Aic[WorstIdx]
      Score <- a
    }
    else if(sortby=="BIC"){
      WorstIdx <- which.max(Bic)
      WorstScore <- Bic[WorstIdx]
      Score <- b
    }
    else{
      WorstIdx <- which.min(LMargLike)
      WorstScore <- - LogMargLike[WorstIdx]
      Score <- - lml
    }

    if(lml > maxlml){
      if(i == 1){
        maxlml <- lml
        totalprob <- 1
      }
      else{
        totalprob <- totalprob * exp(maxlml - lml) + 1
        maxlml <- lml
      }
    }
    else{
      totalprob <- totalprob + exp(lml - maxlml)
    }

    
    if(Score < WorstScore){
      ID[WorstIdx] <- i
      BinId[WorstIdx] <- paste(binid, collapse="")
      Aic[WorstIdx] <- a
      Bic[WorstIdx] <- b
      LogMargLike[WorstIdx] <- lml
      Vars[WorstIdx] <- paste(colnames(Xm),collapse=" ")
    }
    
    ##save full model info, if wanted
    if(storemodels == TRUE){
      gout$BinID <- binid
      gout$ID <- i 
      models[[WorstIdx]] <- gout ##full model information for model i
    }
    if(printi==TRUE)
      print(i)
    i <-  i + 1
  }

  Ar <- rank(Aic) ##create rankings by AIC
  Br <- rank(Bic) ##create rankings by BIC
  Mr <- rank(-LogMargLike) ##create rankings by Marginal Likelihood

  Prob <- exp(LogMargLike - maxlml) / totalprob
  OtherProb <-  1 - sum(Prob)
  

  ##create output list, out$list contains data frame of model selection info
  ##out$models contains model specific information - coefficients, resids, et
  out <- list()
  out$OtherProb <- OtherProb
  out$list <- data.frame(ID=ID, BinaryID=BinId, AIC=Aic, AICrank=Ar, BIC=Bic,
                         BICrank=Br, LogMargLike=LogMargLike, Prob=Prob, MLrank=Mr, Variables=Vars)
  if(storemodels == TRUE)
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
##note: only calculates marginal likelihood up to a multiplicative constant
##that is common to all models
##note2: numerical issue: returns Inf far too easily
marglik <- function(n, k, g, Rsq, Ynorm){
  ##if we have the null model, we know Rsq = 0
  if(k==0)
    Rsq <- 0
  out <- (1+g)/(1+ g*(1-Rsq) )
  out <- out^((n-1)/2)
  out <- out / (1+g)^(k/2)
  ## out <- C * out  ##this is where the common multiplicative constant would come in
  return(out)
}

##Calculates the log marginal likelihood for a given model
##Formula from Liang et al 2007, pg 6, assumes constant g
##later it might be worth implementing with priors chosen for g
##note: only calculates log marginal likelihood up to an additive constant
##that is common to all models
logmarglik <- function(n, k, g, Rsq, Ynorm){
  ##if we have the null model, we know Rsq = 0
  if(k==0)
    Rsq <- 0
  out <- log((1+g)/(1+ g*(1-Rsq) ))
  out <- out * ( (n-1)/2 )
  out <- out - (k/2) * log( (1+g) )
  ## out <- C + out  ##this is where the common additive constant would come in
  return(out)
}

##common multiplicative constant on all marginal likelihoods
##note: not used since it's unnecessary and gamma() blows up for large n
margC <- function(n, Ynorm){
  Ynorm <- sqrt( Ynorm  )^(-n+1)
  K <- gamma((n-1)/2) / pi^((n-1)/2) / sqrt(n)
  out <- K*Ynorm
  return(out)
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
