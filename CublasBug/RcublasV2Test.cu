#include <stdio.h>
#include <stdlib.h>
#include <R.h>
#include <cublas_v2.h>

extern "C" void RcublasV2(int *n){
  float * dX;
  cublasHandle_t han;

  cublasCreate(&han);

  cudaMalloc((void **)&dX, (*n)*sizeof(float));

  cudaFree(dX);

  cublasDestroy(han);


}
