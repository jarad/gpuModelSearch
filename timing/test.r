source("timingfun.r")

n <- 100000
k <- 10


set.seed(2.42)
X <- genX(n,k)
Y <- genY(n,k,X)


time <- system.time(test <- LMsearch(X,Y))

time
