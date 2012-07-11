##error from R:
##Error in gpuLm.fit(x, y, NULL, offset = offset, useSingle = useSingle,  :
##cublas error : getQRDecompBlocked, postblock: : GPU program failed to execute


source("~/gpuModelSearch/Cwrapper/clmsearch.r")
library(gputools)
set.seed(1214)
fittime <- data.frame(n=0, k=0, usrtime=0, systime=0, elapstime=0)

ns <- c(100, 1000, 5000, 10000, 100000, 1000000)
ks <- c(10, 11, 12, 13, 5)

n <- ns[1]
k <- ks[1]

##initial fit because the first time always takes longer than normal for some reason
X <- genX(n,k)
Y <- genY(n,k,X)
system.time(test <- gpuclmsearch(X,Y))

for(i in 1:10){
  if(i != 1){ ##already generated above
    X <- genX(n,k)
    Y <- genY(n,k,X)
  }
  time <- system.time(test <- gpuclmsearch(X,Y))
  fittime[i,] <- c(n,k,time[1:3])
  print(i)
}
 
write.csv(fittime, "cfittime.csv", row.names=F)
fittime <- read.csv("cfittime.csv")

set.seed(310)
k <- ks[1]
n <- ns[2]
##initial fit because the first time always takes longer than normal for some reason
X <- genX(n,k)
Y <- genY(n,k,X)
system.time(test <- gpuclmsearch(X,Y))

for(i in 11:20){
  if(i != 1){ ##already generated above
    X <- genX(n,k)
    Y <- genY(n,k,X)
  }
  time <- system.time(test <- gpuclmsearch(X,Y))
  fittime[i,] <- c(n,k,time[1:3])
  print(i)
}

write.csv(fittime, "cfittime.csv", row.names=F)

fittime <- read.csv("cfittime.csv")

set.seed(333)
n <- ns[3]
k <- ks[1]
##initial fit because the first time always takes longer than normal for some reason
X <- genX(n,k)
Y <- genY(n,k,X)
system.time(test <- gpuclmsearch(X,Y))

for(i in 21:30){
  if(i != 1){ ##already generated above
    X <- genX(n,k)
    Y <- genY(n,k,X)
  }
  time <- system.time(test <- gpuclmsearch(X,Y))
  fittime[i,] <- c(n,k,time[1:3])
  print(i)
}


write.csv(fittime, "cfittime.csv", row.names=F)

fittime <- read.csv("cfittime.csv")

set.seed(345)
n <- ns[4]
k <- ks[1]
##initial fit because the first time always takes longer than normal for some reason
X <- genX(n,k)
Y <- genY(n,k,X)
system.time(test <- gpuclmsearch(X,Y))

for(i in 31:40){
  if(i != 1){ ##already generated above
    X <- genX(n,k)
    Y <- genY(n,k,X)
  }
  time <- system.time(test <- gpuclmsearch(X,Y))
  fittime[i,] <- c(n,k,time[1:3])
  print(i)
}

write.csv(fittime, "cfittime.csv", row.names=F)
fittime <- read.csv("cfittime.csv")

set.seed(351)
n <- ns[5]
k <- ks[1]
##initial fit because the first time always takes longer than normal for some reason
X <- genX(n,k)
Y <- genY(n,k,X)
system.time(test <- gpuLm.fit(X,Y))

for(i in 41:50){
  if(i != 1){ ##already generated above
    X <- genX(n,k)
    Y <- genY(n,k,X)
  }
  time <- system.time(test <- gpuclmsearch(X,Y))
  fittime[i,] <- c(n,k,time[1:3])
  print(i)
}

write.csv(fittime, "cfittime.csv", row.names=F)
fittime <- read.csv("cfittime.csv")

set.seed(424)
n <- ns[6]
k <- ks[1]
##initial fit because the first time always takes longer than normal for some reason
X <- genX(n,k)
Y <- genY(n,k,X)
system.time(test <- gpuLm.fit(X,Y))

for(i in 51:60){
  if(i != 1){ ##already generated above
    X <- genX(n,k)
    Y <- genY(n,k,X)
  }
  time <- system.time(test <- gpuclmsearch(X,Y))
  fittime[i,] <- c(n,k,time[1:3])
  print(i)
}

write.csv(fittime, "cfittime.csv", row.names=F)

fittime <- read.csv("cfittime.csv")


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
  if(i>66){
    time <- system.time(test <- gpuclmsearch(X,Y))
    fittime[i,] <- c(n,k,time[1:3])
  }
  print(i)
}


write.csv(fittime, "cfittime.csv", row.names=F)

fittime <- read.csv("cfittime.csv")

n <- ns[1]
k <- ks[3]
set.seed(358)
##initial fit because the first time always takes longer than normal for some reason
X <- genX(n,k)
Y <- genY(n,k,X)
system.time(test <- gpuLm.fit(X,Y))

for(i in 71:80){
  if(i != 1){ ##already generated above
    X <- genX(n,k)
    Y <- genY(n,k,X)
  }
  if(i>79){
    time <- system.time(test <- gpuclmsearch(X,Y))
    fittime[i,] <- c(n,k,time[1:3])
  }
  print(i)
}

write.csv(fittime, "cfittime.csv", row.names=F)

fittime <- read.csv("cfittime.csv")
n <- ns[1]
k <- ks[4]
set.seed(410)
##initial fit because the first time always takes longer than normal for some reason
X <- genX(n,k)
Y <- genY(n,k,X)
system.time(test <- gpuLm.fit(X,Y))

for(i in 81:90){
  if(i != 1){ ##already generated above
    X <- genX(n,k)
    Y <- genY(n,k,X)
  }
  if(i==90){
    time <- system.time(test <- gpuclmsearch(X,Y))
    fittime[i,] <- c(n,k,time[1:3])
  }
  print(i)
}

write.csv(fittime, "cfittime.csv", row.names=F)

fittime <- read.csv("cfittime.csv")
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
  time <- system.time(test <- gpuclmsearch(X,Y))
  fittime[i,] <- c(n,k,time[1:3])
  print(i)
}

write.csv(fittime, "cfittime.csv", row.names=F)

fittime <- read.csv("cfittime.csv")
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
  time <- system.time(test <- gpuclmsearch(X,Y))
  fittime[i,] <- c(n,k,time[1:3])
  print(i)
}

write.csv(fittime, "cfittime.csv", row.names=F)

fittime <- read.csv("cfittime.csv")
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
  time <- system.time(test <- gpuclmsearch(X,Y))
  fittime[i,] <- c(n,k,time[1:3])
  print(i)
}

write.csv(fittime, "cfittime.csv", row.names=F)

fittime <- read.csv("cfittime.csv")
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
  time <- system.time(test <- gpuclmsearch(X,Y))
  fittime[i,] <- c(n,k,time[1:3])
  print(i)
}

write.csv(fittime, "cfittime.csv", row.names=F)

fittime <- read.csv("cfittime.csv")
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
  time <- system.time(test <- gpuclmsearch(X,Y))
  fittime[i,] <- c(n,k,time[1:3])
  print(i)
}

write.csv(fittime, "cfittime.csv", row.names=F)

fittime <- read.csv("cfittime.csv")
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
  time <- system.time(test <- gpuclmsearch(X,Y))
  fittime[i,] <- c(n,k,time[1:3])
  print(i)
}

write.csv(fittime, "cfittime.csv", row.names=F)
