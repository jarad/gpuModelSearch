#include <stdio.h>
#include <stdlib.h>
#include <cublas.h>
#include <R.h>

extern "C" void Rcublas(int *n){
  float * dX;

  cublasInit();

  cublasAlloc( (*n), sizeof(float), (void **)&dX);

  //checkCublasError("Cublas: ");
  //checkCudaError("Cuda: ");

  cublasFree(dX);

  

  cublasShutdown();

}
