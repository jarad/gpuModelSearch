source("~/gpuModelSearch/Rwrapper/lmsearch.r")

set.seed(1214)
fittime <- data.frame(n=0, k=0, usrtime=0, systime=0, elapstime=0)

n <- 100
k <- 10

for(i in 1:10){
  print(i)
  X <- genX(n,k)
  Y <- genY(n,k,X)
  time <- system.time(test <- gpuLmsearch(Y,X))
  fittime[i,] <- c(n,k,time[1:3])
}

X <- matrix(c(rep(1,10), 1:10), ncol=2)
colnames(X) <- c("Int", "x")
Y <- X%*%t(t(c(1,1)))

test <- gpuLm.fit(X,Y)

dat <- data.frame(X,Y)

test <- gpuLm(data=dat, Y~x)

##error from R:
##Error in gpuLm.fit(x, y, NULL, offset = offset, useSingle = useSingle,  :
##cublas error : getQRDecompBlocked, postblock: : GPU program failed to execute
