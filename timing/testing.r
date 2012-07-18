source("~/gpuModelSearch/Rwrapper/rlmsearch.r")
source("~/gpuModelSearch/Cwrapper/clmsearch.r")
source("~/gpuModelSearch/Csmart/cslmsearch.r")

set.seed(1234)

##to initialize the gpu since this takes ~10 seconds the first time
system.time(test <- gpuLm.fit(matrix(c(1,2,6,9),2,2),matrix(c(1,3),ncol=1)))

n <- 10
k <- 2
X <- genX(n,k)
Y <- genY(n,k,X)

rtime <- system.time(rlm <- gpuRlmsearch(X,Y))
ctime <- system.time(clm <- gpuClmsearch(X,Y))
cstime <- system.time(cslm <- gpuCSlmsearch(X,Y))


