##error from R:
##Error in gpuLm.fit(x, y, NULL, offset = offset, useSingle = useSingle,  :
##cublas error : getQRDecompBlocked, postblock: : GPU program failed to execute

fittime <- data.frame(n=0, k=0, usrtime=0, systime=0, elapstime=0)

source("~/gpuModelSearch/Rwrapper/rlmsearch.r")
set.seed(1214)
fittime <- read.csv("rfittime.csv")
ns <- c(100, 1000, 5000, 10000, 100000, 250000, 500000, 1000000, 1250000, 1500000, 2000000)
ks <- c(10, 11, 12, 13, 5)
n <- ns[1]
k <- ks[1]
##initial fit because the first time always takes longer than normal for some reason
system.time(test <- gpuLm.fit(matrix(c(1,2,6,9),2,2),matrix(c(1,3),ncol=1)))

for(i in 1:5){
  X <- genX(n,k)
  Y <- genY(n,k,X)
  if(i == 5){
    time <- system.time(test <- gpurlmsearch(X,Y))
    fittime[i,] <- c(n,k,time[1:3])
  }
  print(i)
}
write.csv(fittime, "rfittime.csv", row.names=F)

source("~/gpuModelSearch/Rwrapper/rlmsearch.r")
fittime <- read.csv("rfittime.csv")
ns <- c(100, 1000, 5000, 10000, 100000, 250000, 500000, 1000000, 1250000, 1500000, 2000000)
ks <- c(10, 11, 12, 13, 5)
set.seed(310)
k <- ks[1]
n <- ns[2]
##initial fit because the first time always takes longer than normal for some reason
system.time(test <- gpuLm.fit(matrix(c(1,2,6,9),2,2),matrix(c(1,3),ncol=1)))

for(i in 6:10){
  X <- genX(n,k)
  Y <- genY(n,k,X)
  if(i == 10){
    time <- system.time(test <- gpurlmsearch(X,Y))
    fittime[i,] <- c(n,k,time[1:3])
  }
  print(i)
}
write.csv(fittime, "rfittime.csv", row.names=F)

source("~/gpuModelSearch/Rwrapper/rlmsearch.r")
fittime <- read.csv("rfittime.csv")
ns <- c(100, 1000, 5000, 10000, 100000, 250000, 500000, 1000000, 1250000, 1500000, 2000000)
ks <- c(10, 11, 12, 13, 5)
set.seed(333)
k <- ks[1]
n <- ns[3]
##initial fit because the first time always takes longer than normal for some reason
system.time(test <- gpuLm.fit(matrix(c(1,2,6,9),2,2),matrix(c(1,3),ncol=1)))

for(i in 11:15){
  X <- genX(n,k)
  Y <- genY(n,k,X)
  if(i==15){
    time <- system.time(test <- gpurlmsearch(X,Y))
    fittime[i,] <- c(n,k,time[1:3])
  }
  print(i)
}
write.csv(fittime, "rfittime.csv", row.names=F)

source("~/gpuModelSearch/Rwrapper/rlmsearch.r")
fittime <- read.csv("rfittime.csv")
ns <- c(100, 1000, 5000, 10000, 100000, 250000, 500000, 1000000, 1250000, 1500000, 2000000)
ks <- c(10, 11, 12, 13, 5)
set.seed(345)
k <- ks[1]
n <- ns[4]
##initial fit because the first time always takes longer than normal for some reason
system.time(test <- gpuLm.fit(matrix(c(1,2,6,9),2,2),matrix(c(1,3),ncol=1)))

for(i in 16:20){
  X <- genX(n,k)
  Y <- genY(n,k,X)
  if(i==20){
    time <- system.time(test <- gpurlmsearch(X,Y))
    fittime[i,] <- c(n,k,time[1:3])
  }
  print(i)
}
write.csv(fittime, "rfittime.csv", row.names=F)

source("~/gpuModelSearch/Rwrapper/rlmsearch.r")
fittime <- read.csv("rfittime.csv")
ns <- c(100, 1000, 5000, 10000, 100000, 250000, 500000, 1000000, 1250000, 1500000, 2000000)
ks <- c(10, 11, 12, 13, 5)
set.seed(351)
k <- ks[1]
n <- ns[5]
##initial fit because the first time always takes longer than normal for some reason
system.time(test <- gpuLm.fit(matrix(c(1,2,6,9),2,2),matrix(c(1,3),ncol=1)))

for(i in 21:25){
  X <- genX(n,k)
  Y <- genY(n,k,X)
  if(i==25){
    time <- system.time(test <- gpurlmsearch(X,Y))
    fittime[i,] <- c(n,k,time[1:3])
  }
  print(i)
}
write.csv(fittime, "rfittime.csv", row.names=F)

source("~/gpuModelSearch/Rwrapper/rlmsearch.r")
fittime <- read.csv("rfittime.csv")
ns <- c(100, 1000, 5000, 10000, 100000, 250000, 500000, 1000000, 1250000, 1500000, 2000000)
ks <- c(10, 11, 12, 13, 5)
set.seed(424)
n <- ns[6]
k <- ks[1]
##initial fit because the first time always takes longer than normal for some reason
system.time(test <- gpuLm.fit(matrix(c(1,2,6,9),2,2),matrix(c(1,3),ncol=1)))

for(i in 26:30){
  X <- genX(n,k)
  Y <- genY(n,k,X)
  if(i==30){
    time <- system.time(test <- gpurlmsearch(X,Y))
    fittime[i,] <- c(n,k,time[1:3])
  }
  print(i)
}
write.csv(fittime, "rfittime.csv", row.names=F)

source("~/gpuModelSearch/Rwrapper/rlmsearch.r")
fittime <- read.csv("rfittime.csv")
ns <- c(100, 1000, 5000, 10000, 100000, 250000, 500000, 1000000, 1250000, 1500000, 2000000)
ks <- c(10, 11, 12, 13, 5)
n <- ns[7]
k <- ks[1]
set.seed(522)
##initial fit because the first time always takes longer than normal for some reason
system.time(test <- gpuLm.fit(matrix(c(1,2,6,9),2,2),matrix(c(1,3),ncol=1)))

for(i in 31:35){
  X <- genX(n,k)
  Y <- genY(n,k,X)
  if(i==35){
    time <- system.time(test <- gpurlmsearch(X,Y))
    fittime[i,] <- c(n,k,time[1:3])
  }
  print(i)
}
write.csv(fittime, "rfittime.csv", row.names=F)


source("~/gpuModelSearch/Rwrapper/rlmsearch.r")
fittime <- read.csv("rfittime.csv")
ns <- c(100, 1000, 5000, 10000, 100000, 250000, 500000, 1000000, 1250000, 1500000, 2000000)
ks <- c(10, 11, 12, 13, 5)
n <- ns[8]
k <- ks[1]
set.seed(358)
##initial fit because the first time always takes longer than normal for some reason
system.time(test <- gpuLm.fit(matrix(c(1,2,6,9),2,2),matrix(c(1,3),ncol=1)))

for(i in 36:40){
  X <- genX(n,k)
  Y <- genY(n,k,X)  
  if(i==40){
    time <- system.time(test <- gpurlmsearch(X,Y))
    fittime[i,] <- c(n,k,time[1:3])
  }
  print(i)
}
write.csv(fittime, "rfittime.csv", row.names=F)

source("~/gpuModelSearch/Rwrapper/rlmsearch.r")
fittime <- read.csv("rfittime.csv")
ns <- c(100, 1000, 5000, 10000, 100000, 250000, 500000, 1000000, 1250000, 1500000, 2000000)
ks <- c(10, 11, 12, 13, 5)
n <- ns[9]
k <- ks[1]
set.seed(410)
##initial fit because the first time always takes longer than normal for some reason
system.time(test <- gpuLm.fit(matrix(c(1,2,6,9),2,2),matrix(c(1,3),ncol=1)))

for(i in 41:45){
  X <- genX(n,k)
  Y <- genY(n,k,X)
  if(i==45){
    time <- system.time(test <- gpurlmsearch(X,Y))
    fittime[i,] <- c(n,k,time[1:3])
  }
  print(i)
}
write.csv(fittime, "rfittime.csv", row.names=F)

fittime <- read.csv("rfittime.csv")
ns <- c(100, 1000, 5000, 10000, 100000, 1000000)
ks <- c(10, 11, 12, 13, 5)
n <- ns[1]
k <- ks[5]
set.seed(446)
##initial fit because the first time always takes longer than normal for some reason
X <- genX(n,k)
Y <- genY(n,k,X)
system.time(test <- gpuLm.fit(X,Y))
for(i in 91:100){
  if(i != 1){ ##already generated above
    X <- genX(n,k)
    Y <- genY(n,k,X)
  }
  time <- system.time(test <- gpurlmsearch(X,Y))
  fittime[i,] <- c(n,k,time[1:3])
  print(i)
}

write.csv(fittime, "rfittime.csv", row.names=F)

fittime <- read.csv("rfittime.csv")
ns <- c(100, 1000, 5000, 10000, 100000, 1000000)
ks <- c(10, 11, 12, 13, 5)
n <- ns[2]
k <- ks[5]
set.seed(448)
##initial fit because the first time always takes longer than normal for some reason
X <- genX(n,k)
Y <- genY(n,k,X)
system.time(test <- gpuLm.fit(X,Y))
for(i in 101:110){
  if(i != 1){ ##already generated above
    X <- genX(n,k)
    Y <- genY(n,k,X)
  }
  time <- system.time(test <- gpurlmsearch(X,Y))
  fittime[i,] <- c(n,k,time[1:3])
  print(i)
}

write.csv(fittime, "rfittime.csv", row.names=F)

fittime <- read.csv("rfittime.csv")
ns <- c(100, 1000, 5000, 10000, 100000, 1000000)
ks <- c(10, 11, 12, 13, 5)
n <- ns[3]
k <- ks[5]
set.seed(449)
##initial fit because the first time always takes longer than normal for some reason
X <- genX(n,k)
Y <- genY(n,k,X)
system.time(test <- gpuLm.fit(X,Y))
for(i in 111:120){
  if(i != 1){ ##already generated above
    X <- genX(n,k)
    Y <- genY(n,k,X)
  }
  time <- system.time(test <- gpurlmsearch(X,Y))
  fittime[i,] <- c(n,k,time[1:3])
  print(i)
}

write.csv(fittime, "rfittime.csv", row.names=F)

fittime <- read.csv("rfittime.csv")
ns <- c(100, 1000, 5000, 10000, 100000, 1000000)
ks <- c(10, 11, 12, 13, 5)
n <- ns[4]
k <- ks[5]
set.seed(450)
##initial fit because the first time always takes longer than normal for some reason
X <- genX(n,k)
Y <- genY(n,k,X)
system.time(test <- gpuLm.fit(X,Y))
for(i in 121:130){
  if(i != 1){ ##already generated above
    X <- genX(n,k)
    Y <- genY(n,k,X)
  }
  time <- system.time(test <- gpurlmsearch(X,Y))
  fittime[i,] <- c(n,k,time[1:3])
  print(i)
}

write.csv(fittime, "rfittime.csv", row.names=F)

fittime <- read.csv("rfittime.csv")
ns <- c(100, 1000, 5000, 10000, 100000, 1000000)
ks <- c(10, 11, 12, 13, 5)
n <- ns[5]
k <- ks[5]
set.seed(451)
##initial fit because the first time always takes longer than normal for some reason
X <- genX(n,k)
Y <- genY(n,k,X)
system.time(test <- gpuLm.fit(X,Y))
for(i in 131:140){
  if(i != 1){ ##already generated above
    X <- genX(n,k)
    Y <- genY(n,k,X)
  }
  time <- system.time(test <- gpurlmsearch(X,Y))
  fittime[i,] <- c(n,k,time[1:3])
  print(i)
}

write.csv(fittime, "rfittime.csv", row.names=F)

fittime <- read.csv("rfittime.csv")
ns <- c(100, 1000, 5000, 10000, 100000, 1000000)
ks <- c(10, 11, 12, 13, 5)
n <- ns[6]
k <- ks[5]
set.seed(453)
##initial fit because the first time always takes longer than normal for some reason
X <- genX(n,k)
Y <- genY(n,k,X)
system.time(test <- gpuLm.fit(X,Y))
for(i in 141:150){
  if(i != 1){ ##already generated above
    X <- genX(n,k)
    Y <- genY(n,k,X)
  }
  time <- system.time(test <- gpurlmsearch(X,Y))
  fittime[i,] <- c(n,k,time[1:3])
  print(i)
}

write.csv(fittime, "rfittime.csv", row.names=F)
