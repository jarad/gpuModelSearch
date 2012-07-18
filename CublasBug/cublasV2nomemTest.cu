#include <stdio.h>
#include <stdlib.h>
#include <cublas_v2.h>

int main(void){
  int i;
  for(i = 0; i < 100000; i++){
    cublasHandle_t han;
    
    printf(" %d \n", i);

    cublasCreate(&han);
    
    cublasDestroy(han);

  }
  return 0;
}
