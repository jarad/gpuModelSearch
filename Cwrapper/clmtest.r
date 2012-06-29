dyn.load("lmsearch.so")
X <- matrix(c(rep(1,5),1:5,(1:5)^2),ncol=3)
Y <- rnorm(5)+X%*%t(t(c(1,2,3)))
g <- 5L
num_save <- 3L
aics <- rep(10000000000, num_save)
mode(aics) <- "single"
bics <- aics
logmargs <- -aics/100
models <- integer(num_save)
mode(Y) <- "single"
mode(X) <- "single"
sorttype <- 2L

out <- .C("lmsearch",X,as.integer(5),as.integer(3),Y,as.integer(1),g,aic=aics,bic=bics,lml=logmargs,id=models,num_save,sorttype)

source("~/gpuModelSearch/Rwrapper/rlmsearch.r")

out2 <- gpuLmsearch(X, Y, sortby="BIC")
