##This script reproduces a bug in gputools. After calling gpuLm.fit enough times in
##one session, I get the following error on any call to gpuLm.fit:

##Error in gpuLm.fit(X, Y) :
##  cublas error : getQRDecompBlocked, postblock: : GPU program failed to execute

##Closing the R process and reopening it solves the problem. The following script
##reproduces the error as well as producing another potential bug that is probably
##related. We think the issue is a memory leak somewhere.

library(gputools)

n <- 1000 ##sample size
k <- 200  ##number of covariates

##create model matrix and response vector
set.seed(1)
X <- matrix(c(rep(1,n), rnorm(k*n)), ncol=k+1)
Y <- X %*% t(t( c( 1:(k+1) ) )) + rnorm(n)

##for storing how long it takes to fit the model - this will be interesting
fittime <- data.frame(i=0, time=0)

##Fit the model over and over again and store how long it takes to fit each time.
##After fitting the model enough times, R should throw the error above.
##R is still usable at this point, but calls to the GPU are not.
M <- 15000
for(i in 1:M){
  fittime[i,] <- c(i, system.time(gpuLm.fit(X,Y))[3])
  if(i %% 100 == 0){
    print(i)
  }
}

##You may need to increase M to reproduce the error, it occurs for me at i=13173
##Increasing n or k seems to decrease the number of iterations before the error.

##now from my previous work, it appeared that as i increased, the amount of time
##it took to fit a model also increased. This is much more visible when I'm timing
##how long it takes to fit all possible sub-models, but linear regression picks it
##up in this case:

o <- lm(data=fittime[-1,], time~i)
summary(o)

##I remove the first observation because the first call to the gpu from a given process
##always takes longer for some reason on our system - whether it's R or just a
##cuda C program.

##here's the table of coefficients I get:
##Coefficients:
##             Estimate Std. Error t value Pr(>|t|)
##(Intercept) 1.401e-01  9.112e-05  1538.1   <2e-16 ***
##i           2.360e-07  1.198e-08    19.7   <2e-16 ***

##The slope is really small but positive and highly significant. So each time we fit
##the model, it takes slightly longer to fit. The increase is visible to the naked
##eye if I time how long it takes to fit all possible submodels. I bet it would also
##be visible for large enough n and k.

##So two potential bugs that I think are related:
##1) After calling gpuLm.fit enough times, the gpu will no longer respond.
##2) If we fit the same model using gpuLm.fit over and over again in the same R
## process, it tends to take longer to fit the model as the number of iterations
## increases.
