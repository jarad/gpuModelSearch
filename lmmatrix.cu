/*
  This piece of code takes a full model matrix X with p columns (covariates) and
  n rows (observations) and computes the model matrix for each possible 
  submodel including at least 1 column (covariate) from X.

  issues:
    - Currently each model, indexed by id, gets a matrix of size n*p where the 
      cols corresponding to variables NOT in model id are set to 0. Instead, each
      model should get the correct n*ncol[id] matrix with those columns omitted.
    - Currently the size of matrix X must be set at compile time. Dynamic allocation
      needs to be implemented so that this can be set an run time.
 */


#include <stdio.h>
#include <stdlib.h>
#include <math.h>


//Converts the model ID into a binary string of variable (column) indicators
__device__ void modelid(int id, int p, int *var, int * ncol){
  
  int i, k; //loop counters
  int remain = id+1; //remainder - it's + 1 because model 0 has no model matrix
  int tmpcol = 0; //column counter

  //create the binary variable indicator and count number of columns in model matrix
  for(i=0; i<p; i++){
    k = 1 << (p - i - 1);
    if (remain >= k){
      var[i+id*p] = 1;
      remain = remain % k;
      tmpcol += 1;
    }
    else{
      var[i+id*p] = 0;
    }
  }

  ncol[id] = tmpcol;

}


//Creates the model matrix for each possible model
__global__ void modelmatrix(int * p, int * n, int * N, float * X, float * Xm, 
			    int * ncol, int * var){
  int id = blockIdx.x; //model id
  int i, j; //loop counters
  
  if(id < *N){
    //convert model ID (id) to binary variable list
    modelid(id, *p, var, ncol);
    
    //initialize model id's model matrix
    for(i=0; i<(*n); i++)
      for(j=0; j< (*p); j++)
	if( var[j + id*(*p)] == 1 )
	  Xm[ id + (i*(*p)+j)*(*N) ] = X[j + i*(*p)];
	else
	  Xm[ id + (i*(*p)+j)*(*N) ] = 0;
  }
}



int main (void){
  int p = 3;  //# covariates including intercept, i.e. cols of full model matrix
  int N = (1 << p) - 1; //Number of possible models
  int n=5;  //# observations, i.e. rows of model matrix
  int * dp, * dn, * dN;  //device p and device n
  int i, j, k;  //loop counters
  float X[n][p]; //full model matrix
  float Xm[n][p][N]; //array of all possible model matricies
  float * dX, * dXm; //device X and device Xm
  int ncol[N]; //# columns in each model's model matrix
  int * dncol; //device ncol
  int var[N][p]; //binary string indicating which variables are in which models
  int * dvar;  //device var

  //initialize X
  for(i=0; i<n; i++)
    for(j=0; j<p; j++)
      X[i][j]=pow(i+1,j);

  //allocate memory on device for device variables
  cudaMalloc( (void**)&dp, sizeof(int));
  cudaMalloc( (void**)&dn, sizeof(int));
  cudaMalloc( (void**)&dN, sizeof(int));
  cudaMalloc( (void**)&dX, p*n*sizeof(float));
  cudaMalloc( (void**)&dncol, N*sizeof(int));
  cudaMalloc( (void**)&dvar, N*p*sizeof(int));
  cudaMalloc( (void**)&dXm, N*n*p*sizeof(float));
  
  //copy device variables from host
  cudaMemcpy( dp, &p, sizeof(int), cudaMemcpyHostToDevice);
  cudaMemcpy( dn, &n, sizeof(int), cudaMemcpyHostToDevice);
  cudaMemcpy( dN, &N, sizeof(int), cudaMemcpyHostToDevice);
  cudaMemcpy( dX, X, p*n*sizeof(float), cudaMemcpyHostToDevice);

  //Create each model's model matrix
  modelmatrix<<<N,1>>>(dp, dn, dN, dX, dXm, dncol, dvar);

  //copy updated device variables to host
  cudaMemcpy(ncol, dncol, N*sizeof(int), cudaMemcpyDeviceToHost);
  cudaMemcpy(var, dvar, p*N*sizeof(int), cudaMemcpyDeviceToHost);
  cudaMemcpy(Xm, dXm, p*N*n*sizeof(float),cudaMemcpyDeviceToHost);

  //free memory
  cudaFree(dvar);
  cudaFree(dncol);
  cudaFree(dp);
  cudaFree(dn);
  cudaFree(dX);
  cudaFree(dN);

  //Print each model's model matrix
  for(k=0; k<N; k++){
    printf("Model: %d Cols: %d Var: ", k+1, ncol[k]);
    for(i=0; i<p; i++){
      printf("%d",var[k][i]);
    }
    printf("\n\nModelMatrix:\n\n");
    
    for(i=0; i<n; i++){
	for(j=0; j<p; j++){
	  printf("%.0f  ", Xm[i][j][k] );
	}
	printf("\n");
      }
    printf("\n\n");
  }

  //print the full model matrix (for debugging)
  //  printf("Full Matrix X\n\n");
  //  for(i=0; i<n; i++){
  //    for(j=0; j<p; j++){
  //      printf("%.f ", X[i][j]);
  //    }
  //    printf("\n");
  //  }
	  
  return 0;
}
