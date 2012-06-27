##error from R:
##Error in gpuLm.fit(x, y, NULL, offset = offset, useSingle = useSingle,  :
##cublas error : getQRDecompBlocked, postblock: : GPU program failed to execute


source("~/gpuModelSearch/Rwrapper/rlmsearch.r")

set.seed(1214)
fittime <- data.frame(n=0, k=0, usrtime=0, systime=0, elapstime=0)

ns <- c(100, 1000, 5000, 10000, 100000, 1000000)
ks <- c(10, 15, 20, 50, 100, 500, 1000)

n <- ns[1]
k <- ks[1]

##initial fit because the first time always takes longer than normal for some reason
X <- genX(n,k)
Y <- genY(n,k,X)
system.time(test <- gpuLmsearch(X,Y))

for(i in 1:10){
  if(i != 1){ ##already generated above
    X <- genX(n,k)
    Y <- genY(n,k,X)
  }
  time <- system.time(test <- gpuLmsearch(X,Y))
  fittime[i,] <- c(n,k,time[1:3])
  print(i)
}

write.csv(fittime, "fittime.csv", row.names=F)
fittime <- read.csv("fittime.csv")

set.seed(310)
n <- ns[2]
##initial fit because the first time always takes longer than normal for some reason
X <- genX(n,k)
Y <- genY(n,k,X)
system.time(test <- gpuLmsearch(X,Y, printi=TRUE))

for(i in 11:20){
  if(i != 1){ ##already generated above
    X <- genX(n,k)
    Y <- genY(n,k,X)
  }
  time <- system.time(test <- gpuLmsearch(X,Y))
  fittime[i,] <- c(n,k,time[1:3])
  print(i)
}

write.csv(fittime, "fittime.csv", row.names=F)
fittime <- read.csv("fittime.csv")

set.seed(333)
n <- ns[3]
##initial fit because the first time always takes longer than normal for some reason
X <- genX(n,k)
Y <- genY(n,k,X)
system.time(test <- gpuLmsearch(X,Y, printi=TRUE))

for(i in 21:30){
  if(i != 1){ ##already generated above
    X <- genX(n,k)
    Y <- genY(n,k,X)
  }
  time <- system.time(test <- gpuLmsearch(X,Y))
  fittime[i,] <- c(n,k,time[1:3])
  print(i)
}



write.csv(fittime, "fittime.csv", row.names=F)
fittime <- read.csv("fittime.csv")

set.seed(345)
n <- ns[4]
##initial fit because the first time always takes longer than normal for some reason
X <- genX(n,k)
Y <- genY(n,k,X)
system.time(test <- gpuLmsearch(X,Y, printi=TRUE))

for(i in 31:40){
  if(i != 1){ ##already generated above
    X <- genX(n,k)
    Y <- genY(n,k,X)
  }
  time <- system.time(test <- gpuLmsearch(X,Y))
  fittime[i,] <- c(n,k,time[1:3])
  print(i)
}

write.csv(fittime, "fittime.csv", row.names=F)
fittime <- read.csv("fittime.csv")

set.seed(351)
n <- ns[5]
##initial fit because the first time always takes longer than normal for some reason
X <- genX(n,k)
Y <- genY(n,k,X)
system.time(test <- gpuLm.fit(X,Y))

for(i in 41:50){
  if(i != 1){ ##already generated above
    X <- genX(n,k)
    Y <- genY(n,k,X)
  }
  time <- system.time(test <- gpuLmsearch(X,Y))
  fittime[i,] <- c(n,k,time[1:3])
  print(i)
}

write.csv(fittime, "fittime.csv", row.names=F)
fittime <- read.csv("fittime.csv")

set.seed(424)
n <- ns[6]
##initial fit because the first time always takes longer than normal for some reason
X <- genX(n,k)
Y <- genY(n,k,X)
system.time(test <- gpuLm.fit(X,Y))

for(i in 51:60){
  if(i != 1){ ##already generated above
    X <- genX(n,k)
    Y <- genY(n,k,X)
  }
  time <- system.time(test <- gpuLmsearch(X,Y))
  fittime[i,] <- c(n,k,time[1:3])
  print(i)
}

write.csv(fittime, "fittime.csv", row.names=F)
fittime <- read.csv("fittime.csv")


n <- ns[1]
k <- ks[2]
set.seed(522)
##initial fit because the first time always takes longer than normal for some reason
X <- genX(n,k)
Y <- genY(n,k,X)
system.time(test <- gpuLm.fit(X,Y))

for(i in 61:70){
  if(i != 1){ ##already generated above
    X <- genX(n,k)
    Y <- genY(n,k,X)
  }
  time <- system.time(test <- gpuLmsearch(X,Y))
  fittime[i,] <- c(n,k,time[1:3])
  print(i)
}
write.csv(fittime, "fittime.csv", row.names=F)
