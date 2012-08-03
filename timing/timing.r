source("timingfun.r")

fittime <- data.frame(n=0, k=0, usrtime=0, systime=0, elapstime=0, wrap="C")
levels(fittime$wrap) <- c("C", "R", "CS")

ns <- c(100, 1000, 5000, 10000, 100000, 250000, 500000, 1000000, 1250000, 1500000, 2000000)  ## number of rows/obs
ks <- c(10, 11, 12, 13, 5) ## number of columns
wraps <- c("R", "C", "CS") ## which wrapper to use
N <- 5 ## number of replications of each combination of factors

set.seed(3.57)

## Set gpu to one that no one is using
chooseGpu(3)

i <- 1
for(n in ns){
  for(k in ks){
    for(wrap in wraps){
      for(j in 1:N){

        X <- genX(n,k)
        Y <- genY(n,k,X)

        ## just to initialize gpu for timing purposes
        gpuSolve(diag(2))
      
        if(wrap == "R"){
          time <- system.time(test <- gpuRlmsearch(X,Y))
        }
        else if(wrap == "C"){
          time <- system.time(test <- gpuClmsearch(X,Y))
        }
        else if(wrap == "CS"){
          time <- system.time(test <- gpuCSlmsearch(X,Y))
        }
        
        fittime[i,] <- c(n,k,time[1:3], wrap)
        i <- i + 1

        ## reset the gpu here
        gpuReset()
        
      }
    }
  }
}

write.csv(fittime, "fittime.csv", row.names=F)
