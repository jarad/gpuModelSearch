#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<cublas.h>

__global__ void gpuLMSearchFitF(float * X, int rows, int cols, float * Y, int yCols,
				double tol, float * coeffs, float * resids,
				float * effects, float * aic, float * bic, int mods, 
				int * rank, int *pivot, double * drAux){
  /*The basic algorithm:
    
    1. Use thread and / or block id to determine model ID.
    2. Use model ID to determine which columns of X to use.
    3. Use Householder transformation to find the PQR decomp of Xm - 
       this algorithm is used in gputools.
    4. Use PQR decomp of Xm to find coefficients, residuals, and effects - 
       this is implemented in gputools
    5. Calculate AIC, BIC, and other possible model comparison quantities
    6. Return everything back to the host
  */

  int modelID = blockIdx.y
  

}
