#include <stdio.h>
#include <stdlib.h>
#include <R.h>

extern "C" void Rcuda(int *n){
  float * dX;

  cudaMalloc((void **)&dX, (*n)*sizeof(float));

  cudaFree(dX);


}
