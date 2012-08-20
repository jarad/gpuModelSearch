source("timingfun.r")

fittime <- read.csv("fittime.csv")

levels(fittime$wrap) <- c(levels(fittime$wrap), "LM")

ns <- c(100, 1000, 5000)
ks <- c(16, 18, 20)
N <- 5 ## number of replications of each combination of factors

set.seed(2.42)

## Set gpu to one that no one is using
chooseGpu(3)

i <- dim(fittime)[1]
for(n in ns){
  for(k in ks){
    for(j in 1:N){

      X <- genX(n,k)
      Y <- genY(n,k,X)
      
      ## just to initialize gpu for timing purposes
      gpuSolve(diag(2))
      
      time <- system.time(test <- gpuCSlmsearch(X,Y))
      
          
      fittime[i,] <- c(n,k,time[1:3], "CS")
      write.csv(fittime, "fittime.csv", row.names=F)
      
      print(c(i,n,k,"CS",j))
          
      i <- i + 1
          
      ## reset the gpu here
      gpuReset()
      
    }
  }
}

ns <- c(100, 1000, 5000, 10000, 100000, 250000, 500000, 1000000, 1250000, 1500000, 2000000)
ks <- c(10, 11, 12, 13, 5) ## number of columns

for(n in ns){
  for(k in ks){
    if(k<=10 | n<=100000){
      for(j in 1:N){
        
      X <- genX(n,k)
      Y <- genY(n,k,X)
      
      ## just to initialize gpu for timing purposes
      gpuSolve(diag(2))
      
      time <- system.time(test <- LMsearch(X,Y))



      fittime[i,] <- c(n,k,time[1:3], "LM")
      write.csv(fittime, "fittime.csv", row.names=F)
      
      print(c(i,n,k,"LM",j))
          
      i <- i + 1
          
      ## reset the gpu here
      gpuReset()

      }
    }
  }
}






