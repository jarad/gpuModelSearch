

##Stores model selection info of the best 1000 models according to sortby
LMsearch <- function(X, Y, g=nrow(Y), sortby="AIC"){
  p <- ncol(X) #number of regressors, including intercept
  k <- p - 1 #number of possible covariates
  n <- nrow(X) #sample size
  M <- 2^k #number of possible models (all contain intercept)
  nlist <- min(M, 1000)
  
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

    gout <- lm.fit(Xm, Y) #fit model i on gpu using lm
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

  out$Models <- data.frame(ID=ID, BinaryID=BinId, AIC=Aic, AICrank=Ar, BIC=Bic,
                         BICrank=Br, LogMargLike=LogMargLike, PostProb=Prob,
                         MLrank=Mr, Variables=Vars)
  
  out$OtherProb <- OtherProb

  if(sortby=="AIC")
    out$Models <- out$Models[ order(Ar),]
  else if(sortby=="BIC")
    out$Models <- out$Models[ order(Br),]
  else
    out$Models <- out$Models[ order(Mr),]
  
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
