source("clmsearch.r")

X <- matrix(c(rep(1,100),rnorm(3*100)), ncol=4)

colnames(X) <- c("int", paste("x", 1:3, sep=""))
Y <- rnorm(100) + X%*%t(t(c(1,2,3,4)))

test <- gpulmsearch(X,Y)
