#include <stdio.h>
#include <stdlib.h>
#include <cublas_v2.h>

int main(void){
  int i;
  for(i = 0; i < 100000; i++){
    float * dX;
    cublasHandle_t han;
    
    printf(" %d \n", i);

    cublasCreate(&han);
    
    cudaMalloc((void **)&dX, 10*sizeof(float));
    
    cudaFree(dX);
    
    cublasDestroy(han);

  }
  return 0;
}
