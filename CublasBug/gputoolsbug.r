####This script is outdated but is included hear for informative purposes.
####The bug is actually in the cublas library, not gputools. See README for details.


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

##Fit the model over and over again and store how long it takes to fit each time.
##After fitting the model enough times, R should throw the error above.
##R is still usable at this point, but calls to the GPU are not.
M <- 15000
for(i in 1:M){
  gpuLm.fit(X,Y)
  if(i %% 100 == 0){
    print(i)
  }
}

##You may need to increase M to reproduce the error, it occurs for me at i=13173
##and at i=11963 for another use (with a segmentation fault when exiting R). 
##Increasing n or k seems to decrease the number of iterations before the error.


