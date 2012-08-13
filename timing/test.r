source("timingfun.r")

n <- 500000
k <- 20
chooseGpu(3)

set.seed(2.42)
X <- genX(n,k)
Y <- genY(n,k,X)

gpuSolve(diag(2))


time <- system.time(test <- gpuCSlmsearch(X,Y))

time
