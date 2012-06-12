/*
  This piece of code takes a full model matrix X with p columns (covariates) and
  n rows (observations) and computes the model matrix for each possible 
  submodel including at least 1 column (covariate) from X.

  issues:
    - Currently each model, indexed by id, gets a matrix of size n*p where the 
      cols corresponding to variables NOT in model id are set to 0. Instead, each
      model should get the correct n*ncol[id] matrix with those columns omitted.
    - Storage of each model's matrix (i.e. Xm) is not optimized, but this part will be handled
      differently in the actual fitting algorithm

 */


#include <stdio.h>
#include <stdlib.h>
#include <math.h>


//Converts the model ID into a binary string of variable (column) indicators
__device__ void modelid(int id, int p, int *var, int * ncol, size_t varpitch){
  
  int i, k; //loop counters
  int remain = id+1; //remainder - it's + 1 because model 0 has no cols and is omitted
  int tmpcol = 0; //column counter
  int * varrow = (int *)((char*)var + id*varpitch); //point to model id's binary rep

  //create the binary variable indicator and count number of columns in model matrix
  for(i=0; i<p; i++){
    k = 1 << (p - i - 1);
    if (remain >= k){
      varrow[i] = 1;
      remain = remain % k;
      tmpcol += 1;
    }
    else{
      varrow[i] = 0;
    }
  }

  ncol[id] = tmpcol;

}


//Creates the model matrix for each possible model
__global__ void modelmatrix(int p, int n, int N, float * X, size_t Xpitch,  
			    float * Xm, int * ncol, int * var, size_t varpitch){
  int id = blockIdx.x; //model id
  int i, j; //loop counters
  int * varrow = (int *)((char*)var + id*varpitch); //pointer to model id's binary rep

  if(id < N){
    //convert model ID (id) to binary variable list
    modelid(id, p, var, ncol, varpitch);
    
    //initialize model id's model matrix
    for(i=0; i<n; i++){
      float * Xrow = (float *)((char*)X + i*Xpitch); //pointer to row i of X
      for(j=0; j<p; j++){
    	if( varrow[j] == 1 )
    	  Xm[id + N*(j + i*p)] = Xrow[j];
	else
	  Xm[id + N*(j + i*p)] = 0;
      }
    }
  }
}



int main (void){
  int p = 3;  //# covariates including intercept, i.e. cols of full model matrix
  int N = (1 << p) - 1; //Number of possible models
  int n=5;  //# observations, i.e. rows of model matrix
  int i, j, k;  //loop counters
  float * X; //full model matrix
  float * Xm; //all possible model matricies
  float * dX, * dXm; //device X and device Xm
  int * ncol; //# columns in each model's model matrix
  int * dncol; //device ncol
  int * var; //binary string indicating which variables are in which models
  int * dvar;  //device var
  size_t Xpitch, varpitch; //pitch variables for dX and dvar


  //allocate memory for X, Xm, var, and ncol
  X = (float *) malloc(n*p*sizeof(float));
  Xm = (float *) malloc(N*n*p*sizeof(float));
  var = (int *) malloc(N*p*sizeof(int));
  ncol = (int *) malloc(N*sizeof(int));

  //initialize X
  for(i=0; i<n; i++)
    for(j=0; j<p; j++)
      X[j + i*p]=pow(i+1,j);

  //allocate memory on device for device variables
  //cudaMalloc( (void**)&dX, p*n*sizeof(float));
  cudaMalloc( (void**)&dncol, N*sizeof(int));
  //cudaMalloc( (void**)&dvar, N*p*sizeof(int));
  cudaMalloc( (void**)&dXm, N*n*p*sizeof(float));
  cudaMallocPitch( (void**)&dX, &Xpitch, p*sizeof(float), n);
  cudaMallocPitch( (void**)&dvar, &varpitch, p*sizeof(int), N);
  
  
  //copy device variables from host
  //cudaMemcpy( dX, X, p*n*sizeof(float), cudaMemcpyHostToDevice);
  cudaMemcpy2D(dX, Xpitch, X, p*sizeof(float), p*sizeof(float), n, 
	       cudaMemcpyHostToDevice);

  //Create each model's model matrix
  modelmatrix<<<N,1>>>(p, n, N, dX, Xpitch, dXm, dncol, dvar, varpitch);


  //copy updated device variables to host
  cudaMemcpy(ncol, dncol, N*sizeof(int), cudaMemcpyDeviceToHost);
  //cudaMemcpy(var, dvar, p*N*sizeof(int), cudaMemcpyDeviceToHost);
  cudaMemcpy2D(var, p*sizeof(int), dvar, varpitch, p*sizeof(int), N, 
	       cudaMemcpyDeviceToHost);
  cudaMemcpy(Xm, dXm, p*N*n*sizeof(float), cudaMemcpyDeviceToHost);
  

  //free memory
  cudaFree(dvar);
  cudaFree(dncol);
  cudaFree(dX);
  cudaFree(dXm);

  //Print each model's model matrix
  for(k=0; k<N; k++){
    printf("Model: %d Cols: %d Var: ", k+1, ncol[k]);
    for(i=0; i<p; i++){
      printf("%d",var[i + k*p]);
    }
    printf("\n\nModelMatrix:\n\n");
    
    for(i=0; i<n; i++){
	for(j=0; j<p; j++){
	  printf("%.0f  ", Xm[ k + N*(j + i*p) ] );
	}
	printf("\n");
      }
    printf("\n\n");
  }

  //free memory allocated on host
  free(var);
  free(ncol);
  free(X);
  free(Xm);
			     
			     

  //print the full model matrix (for debugging)
  //  printf("Full Matrix X\n\n");
  //  for(i=0; i<n; i++){
  //    for(j=0; j<p; j++){
  //      printf("%.f ", X[j+p*i]);
  //    }
  //    printf("\n");
  //  }
	  
  return 0;
}
