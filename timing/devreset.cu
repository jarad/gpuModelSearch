#include <stdio.h>
#include <stdlib.h>
#include <R.h>

extern "C" void gpuReset( int * out );

void gpuReset(int * out){
  cudaError_t err;

  err = cudaDeviceReset();

  if(cudaSuccess == err)
    *out = 0;
  else
    *out = 1;
}


