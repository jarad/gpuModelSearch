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

source("~/gpuModelSearch/Rwrapper/rlmsearch.r")
fittime <- read.csv("rfittime.csv")
ns <- c(100, 1000, 5000, 10000, 100000, 250000, 500000, 1000000, 1250000, 1500000, 2000000)
ks <- c(10, 11, 12, 13, 5)
n <- ns[10]
k <- ks[1]
set.seed(446)
##initial fit because the first time always takes longer than normal for some reason
system.time(test <- gpuLm.fit(matrix(c(1,2,6,9),2,2),matrix(c(1,3),ncol=1)))

for(i in 46:50){
  X <- genX(n,k)
  Y <- genY(n,k,X)
  if(i==50){
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
n <- ns[11]
k <- ks[1]
set.seed(448)
##initial fit because the first time always takes longer than normal for some reason
system.time(test <- gpuLm.fit(matrix(c(1,2,6,9),2,2),matrix(c(1,3),ncol=1)))

for(i in 51:55){
  X <- genX(n,k)
  Y <- genY(n,k,X)
  if(i==55){
    time <- system.time(test <- gpurlmsearch(X,Y))
    fittime[i,] <- c(n,k,time[1:3])
  }
  print(i)
}
write.csv(fittime, "rfittime.csv", row.names=F)


###################################

##end of all k=10

###################################


source("~/gpuModelSearch/Rwrapper/rlmsearch.r")
set.seed(12147)
fittime <- read.csv("rfittime.csv")
ns <- c(100, 1000, 5000, 10000, 100000, 250000, 500000, 1000000, 1250000, 1500000, 2000000)
ks <- c(10, 11, 12, 13, 5)
n <- ns[1]
k <- ks[5]
##initial fit because the first time always takes longer than normal for some reason
system.time(test <- gpuLm.fit(matrix(c(1,2,6,9),2,2),matrix(c(1,3),ncol=1)))

for(i in 56:60){
  X <- genX(n,k)
  Y <- genY(n,k,X)
  if(i == 60){
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
set.seed(3107)
k <- ks[5]
n <- ns[2]
##initial fit because the first time always takes longer than normal for some reason
system.time(test <- gpuLm.fit(matrix(c(1,2,6,9),2,2),matrix(c(1,3),ncol=1)))

for(i in 61:65){
  X <- genX(n,k)
  Y <- genY(n,k,X)
  if(i == 65){
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
k <- ks[5]
n <- ns[3]
##initial fit because the first time always takes longer than normal for some reason
system.time(test <- gpuLm.fit(matrix(c(1,2,6,9),2,2),matrix(c(1,3),ncol=1)))

for(i in 66:70){
  X <- genX(n,k)
  Y <- genY(n,k,X)
  if(i==70){
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
set.seed(3457)
k <- ks[5]
n <- ns[4]
##initial fit because the first time always takes longer than normal for some reason
system.time(test <- gpuLm.fit(matrix(c(1,2,6,9),2,2),matrix(c(1,3),ncol=1)))

for(i in 71:75){
  X <- genX(n,k)
  Y <- genY(n,k,X)
  if(i==75){
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
set.seed(3517)
k <- ks[5]
n <- ns[5]
##initial fit because the first time always takes longer than normal for some reason
system.time(test <- gpuLm.fit(matrix(c(1,2,6,9),2,2),matrix(c(1,3),ncol=1)))

for(i in 76:80){
  X <- genX(n,k)
  Y <- genY(n,k,X)
  if(i==80){
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
set.seed(4247)
n <- ns[6]
k <- ks[5]
##initial fit because the first time always takes longer than normal for some reason
system.time(test <- gpuLm.fit(matrix(c(1,2,6,9),2,2),matrix(c(1,3),ncol=1)))

for(i in 81:85){
  X <- genX(n,k)
  Y <- genY(n,k,X)
  if(i==85){
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
k <- ks[5]
set.seed(5227)
##initial fit because the first time always takes longer than normal for some reason
system.time(test <- gpuLm.fit(matrix(c(1,2,6,9),2,2),matrix(c(1,3),ncol=1)))

for(i in 86:90){
  X <- genX(n,k)
  Y <- genY(n,k,X)
  if(i==90){
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
k <- ks[5]
set.seed(3587)
##initial fit because the first time always takes longer than normal for some reason
system.time(test <- gpuLm.fit(matrix(c(1,2,6,9),2,2),matrix(c(1,3),ncol=1)))

for(i in 91:95){
  X <- genX(n,k)
  Y <- genY(n,k,X)  
  if(i==95){
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
k <- ks[5]
set.seed(4107)
##initial fit because the first time always takes longer than normal for some reason
system.time(test <- gpuLm.fit(matrix(c(1,2,6,9),2,2),matrix(c(1,3),ncol=1)))

for(i in 96:100){
  X <- genX(n,k)
  Y <- genY(n,k,X)
  if(i==100){
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
n <- ns[10]
k <- ks[5]
set.seed(4467)
##initial fit because the first time always takes longer than normal for some reason
system.time(test <- gpuLm.fit(matrix(c(1,2,6,9),2,2),matrix(c(1,3),ncol=1)))

for(i in 101:105){
  X <- genX(n,k)
  Y <- genY(n,k,X)
  if(i==105){
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
n <- ns[11]
k <- ks[5]
set.seed(4487)
##initial fit because the first time always takes longer than normal for some reason
system.time(test <- gpuLm.fit(matrix(c(1,2,6,9),2,2),matrix(c(1,3),ncol=1)))

for(i in 106:110){
  X <- genX(n,k)
  Y <- genY(n,k,X)
  if(i==110){
    time <- system.time(test <- gpurlmsearch(X,Y))
    fittime[i,] <- c(n,k,time[1:3])
  }
  print(i)
}
write.csv(fittime, "rfittime.csv", row.names=F)


###################################

##end of all k = 5

###################################


source("~/gpuModelSearch/Rwrapper/rlmsearch.r")
fittime <- read.csv("rfittime.csv")
ns <- c(100, 1000, 5000, 10000, 100000, 250000, 500000, 1000000, 1250000, 1500000, 2000000)
ks <- c(10, 11, 12, 13, 5)
n <- ns[1]
k <- ks[2]
set.seed(448711)
##initial fit because the first time always takes longer than normal for some reason
system.time(test <- gpuLm.fit(matrix(c(1,2,6,9),2,2),matrix(c(1,3),ncol=1)))
for(i in 111:115){
  X <- genX(n,k)
  Y <- genY(n,k,X)
  if(i==115){
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
n <- ns[1]
k <- ks[3]
set.seed(448712)
##initial fit because the first time always takes longer than normal for some reason
system.time(test <- gpuLm.fit(matrix(c(1,2,6,9),2,2),matrix(c(1,3),ncol=1)))
for(i in 116:120){
  X <- genX(n,k)
  Y <- genY(n,k,X)
  if(i==120){
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
n <- ns[1]
k <- ks[4]
set.seed(448713)
##initial fit because the first time always takes longer than normal for some reason
system.time(test <- gpuLm.fit(matrix(c(1,2,6,9),2,2),matrix(c(1,3),ncol=1)))
for(i in 121:125){
  X <- genX(n,k)
  Y <- genY(n,k,X)
  if(i==125){
    time <- system.time(test <- gpurlmsearch(X,Y))
    fittime[i,] <- c(n,k,time[1:3])
  }
  print(i)
}
write.csv(fittime, "rfittime.csv", row.names=F)
