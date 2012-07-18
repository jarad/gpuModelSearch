##See README for instructions on how to use this.


##Rcublas
n <- 10L
i <- 0
if(!is.loaded("Rcublas"))
  dyn.load("~/gpuModelSearch/CublasBug/RcublasTest.so")
while(i < 50000){
  print(i)
  z <- .C("Rcublas", n)
  i <- i + 1
}

##Rcuda
n <- 10L
i <- 0
if(!is.loaded("Rcuda"))
  dyn.load("~/gpuModelSearch/CublasBug/RcudaTest.so")
while(i < 50000){
  print(i)
  z <- .C("Rcuda", n)
  i <- i + 1
}

##RcublasV2
n <- 10L
i <- 0
if(!is.loaded("RcublasV2"))
  dyn.load("~/gpuModelSearch/CublasBug/RcublasV2Test.so")
while(i < 50000){
  print(i)
  z <- .C("RcublasV2", n)
  i <- i + 1
}
