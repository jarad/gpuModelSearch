
library(gputools)

##converts model id to a binary string
modelid <-  function(id, p){
  M <- (2^p) - 1
  binid <- rep(0,p)
  remain <- id
  for(i in 1:p){
    tmp <- 2^(p - i)
    if(remain >= tmp){
      remain <- remain - tmp
      binid[i] <- 1
    }
  }
  return(binid)
}

##takes model id and full model matrix, outputs model matrix for model id
modelmat <- function(X, binid){
  Xm <- matrix(X[,which(binid==1)], ncol=sum(binid))
  colnames(Xm) <- colnames(X)[which(binid==1)]
  return(Xm)
}

##needs to be turned into a function
X <- matrix(c(rep(1,10),rnorm(10)), ncol=2)
colnames(X) <- c("Intercept", "x1")
Y <- rnorm(10) + X%*%t(t(c(2.5, 1.2)))

p <- ncol(X)
n <- nrow(X)
G <- gamma(n/2)
M <- 2^p - 1
models <- list()

ID <- 1:M
BinId <- paste(ID)
Aic <-  ID
Bic <- ID
MargLike <- ID
Vars <- paste(ID)

for(i in ID){
  binid <- modelid(i,p)
  Xm <- modelmat(X, binid)
  pm <- ncol(Xm)
  ##dat <- data.frame(Y,Xm)
  gout <- gpuLm.fit(Xm, Y)
  sighat <- sum(gout$residuals^2)/gout$df.residual
  gout$ID <- i
  gout$BinID <- binid
  models[[i]] <- gout
  BinId[i] <- paste(binid, collapse="")
  Aic[i] <- aic(n,pm,sighat)
  Bic[i] <- bic(n,pm,sighat)
  MargLike[i] <- marglik(n,p,Xm,sighat,G)
  Vars[i] <- paste(colnames(Xm),collapse=" ")
}

out <- list()
out$list <- data.frame(ID=ID, BinaryID=BinId, AIC=Aic, BIC=Bic, MargLike=MargLike,
                       Variables=Vars)
out$models <- models




#note: AIC and BIC defined only up to additive constant
aic <- function(n, p, sighat){
  out <- 2*(p+1) + n*log(sighat)
  return(out)
}

bic <- function(n, p, sighat){
  out <- (p+1)*log(n)+n*log(sighat)
  return(out)
}

##needs to be rewritten for g-prior
marglik <- function(n, p, X, sighat, G){
  ssr <- ( (n-p) * sighat )^(-n/2)
  D <- 1/sqrt( det( t(X)%*%X ) )
  c <- 2^(p/2)*pi^((p-n)/2)
  out <- c*G*ssr*D
  return(out)
}

