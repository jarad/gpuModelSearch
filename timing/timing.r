source("timingfun.r")

fittime <- data.frame(n=0, k=0, usrtime=0, systime=0, elapstime=0, wrap="C"
                      )
ns <- c(100, 1000, 5000, 10000, 100000, 250000, 500000, 1000000, 1250000, 1500000, 2000000)
ks <- c(10, 11, 12, 13, 5)
wraps <- c("R", "C", "CS")

set.seed(3.57)

## rough structure of what needs to happen
i <- 1
for(n in ns){
  for(k in ks){
    for(wrap in wraps){

      X <- genX(n,k)
      Y <- genY(n,k,X)
      
      ## run some simple function that calls the gpu here to initialize it
      
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
      
    }
  }
}

write.csv(fittime, "fittime.csv", row.names=F)
