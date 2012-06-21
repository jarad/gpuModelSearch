#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <cublas.h>
#include <R.h>
#include "cuseful.h"
#include "lsfit.h"
#include "qrdecomp.h"

__global__ void gpuLMSearchFitF(float * X, int rows, int cols, float * Y, int yCols,
				double tol, float * coeffs, float * resids,
				float * effects, float * aic, float * bic, int mods, 
				int * rank, int *pivot, double * drAux){
   

}

float aic(int n, int p, float sighat){
  out = 2*(p+1)+n*log(sighat);
  return(out);
}

float bic(int n, int p, float sighat){
  out = (p+1)*log(n) + n*log(sighat);
  return(out)
}

float LogMargLike(int n, int k, int g, float Rsq, float Ynorm){
  out = log( (1+g) / (1+ g * (1 - Rsq) ) );
  out *= ( (n-1)/2 );
  
}

void lmsearch( float *X, int rows, int cols, float *Y, int ycols, 
			  double tol, double * models ){

  int p = cols, k = cols - 1, M = 1 << k, int n = rows;
  int id, remain, tmp, km, pm, mj;
  int binid[k];
  int i, j;
  float *Xm;
  size_t fbytes = sizeof(float), ibytes = sizeof(int), dbytes = sizeof(double);
  float * coeffs, float *resids, float *effects, int *rank, int *pivot, double *qrAux;
  double tol = 0.0001;

  rank = (int *) malloc(1, ibytes);

  /* For each model:
     1. find the binary ID
     2. convert the binary ID to a matrix
     3. fit the model using the matrix
     4. compute AIC, BIC, and logMargLike
     5. add the model (or not) according to 4
     
     sort output by AIC, BIC, or logMargLike
  */

  for(id=0; id<M; id++){
    km = 0;
    remain = id;
    for(i=0; i<k; i++){
      tmp = 1 << (k - i + 1);
      if( remain >= place ){
	remain -= tmp;
	binid[i] = 1;
	km += 1;
      }
    }
    pm = km + 1;

    Xm = (float *) malloc( n * (km+1) * fbytes );

    mj = 0;

    for(j=0; j<p; j++){    
      if(j == 0 || binid[j-1] == 1){
	for(i=0; i<n; i++){
	  Xm[mj + i*pm] = X[ j + i*p ]
	}
	mj += 1;
      }
    }

    coeffs = (float *) calloc(pm, fbytes);
    resids = (float *) calloc(n, fbytes);

    effects = (float *) malloc(n*ycol*fbytes);
    qrAux = (double *) calloc(pm, dbytes);
    memcpy(effects, Y, n*ycol*fbytes);

    *rank = 1;

    pivot = (int *) malloc(pm, ibytes);
    for(i=0; i<pm; i++)
      pivot[i] = i;

    
    

    gpuLSFitF(Xm, n, pm, Y, ycols, tol, coeffs, resids, effects, rank, pivot, qrAux);

    
    

    free(coeffs);
    free(resids);
    free(effects);
    free(qrAux);
    free(pivot);
    
    free(Xm);
    
  }

  free(rank);

}
