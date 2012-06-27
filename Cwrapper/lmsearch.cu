#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <cublas.h>
#include <R.h>
#include "cuseful.h"
#include "lsfit.h"
#include "qrdecomp.h"

//function prototype for the function to be called from R
extern "C" void lmsearch(float *X, int *rows, int *cols, float *Y, int *ycols, 
			 int *g, float *aics, float *bics, float *logmargs, 
			 int *models, int *num_save, int *sorttype );

//computes aic
float aic(int n, int p, float sighat){
  float out = 2*(p+1)+n*log(sighat);
  return(out);
}

//computes bic
float bic(int n, int p, float sighat){
  float out = (p+1)*log(n) + n*log(sighat);
  return(out);
}

//computes log marginal likelihood
float logmarglike(int n, int k, int g, float Rsq){
  if(k==0){
    float out = 0;
    return(out);
  }
  else{
    float out = (( ((float)n) - 1) / 2) * log((1+g)/(1+ g * (1 - Rsq)));
    out -= ( ((float)k) / 2) * log(1 + g);
    return(out);
  }
}

//finds the worst largest or smallest in src and output the index and the element
void arraysort(float *src, int *idx, float *dest, int num_elem, int max){
  float tmp = *src;
  int id = 0, i;

  for(i = 0; i < num_elem; i++){
    if(max == 1){
      if(tmp < *(src + i) ){
	tmp = *(src + i);
	id = i;
      }
    }
    else{
      if(tmp > *(src + i) ){
	tmp = (*src + i);
	id = i;
      }
    }
    
  }

  *dest = tmp;
  *idx = id;
}

//performs model search
void lmsearch(float *X, int *rows, int *cols, float *Y, int *ycols, int *g, 
	      float *aics, float *bics, float *logmargs, int *models, 
	      int *num_save, int *sorttype ){

  int p = *cols, k = p - 1, M = 1 << k, n = *rows;
  int id, remain, tmp, km, pm, mj;
  int binid[k];  //binary id of current model
  int i, j;
  float *Xm; //pointer to model matrix of current model
  size_t fbytes = sizeof(float), ibytes = sizeof(int), dbytes = sizeof(double);
  float *coeffs, *resids, *effects; //pointers to arrays, needed for gpuLmfit
  int *rank, *pivot; 
  double *qrAux;
  double tol = 0.0001; //tolerance of the fit, float only supported by gputools
  float score, *worstscore; 
  float sighat, ynorm=0, ssr, ymn=0, Rsq;
  int *worstid;
  

  rank = (int *) malloc(ibytes);
  worstid = (int *) malloc(ibytes);
  worstscore = (float *) malloc(fbytes);

  if(rank == NULL){
    Rprintf("\nMemory for rank failed to allocate\n");
    exit(1);
  }

  if(worstid == NULL){
    Rprintf("\nMemory for worstid failed to allocate\n");
    exit(1);
  }

  if(worstscore == NULL){
    Rprintf("\nMemory for worstscore failed to allocate\n");
    exit(1);
  }
  
  //calculate the mean value of the response, Y
  for(i = 0; i < n; i++)
    for(j = 0; j < *ycols; j++)
      ymn += Y[i + j * n];
  ymn /= n;

  //calculate (Y-Ybar)^2, used in marginal likelihood computation
  for(i = 0; i < n; i++)
    for(j = 0; j < *ycols; j++)
      ynorm += pow( Y[i + j * n] - ymn, 2 );

  //for each model
  for(id = 0; id < M; id++){
    //find the binary id and num of covariates in current model
    km = 0; 
    remain = id;
    for(i = 0; i < k; i++){
      tmp = 1 << (k - i - 1);
      if( remain >= tmp ){
	remain -= tmp;
	binid[i] = 1;
	km += 1;
      }
      else{
	binid[i] = 0;
      }
    }
    pm = km + 1;
    
    //allocate memory for model matrix of current model
    Xm = (float *) malloc( n * pm * fbytes );

    if(Xm == NULL){
      Rprintf("\nMemory for Xm failed to allocate\n");
    exit(1);
  }

    mj = 0; //keeps track of column numbers of current model matrix

    //copy the relevant columns of the full model matrix into the current model's
    //model matrix. Note: R stores matrices/arrays in column major format.
    for(j = 0; j < p; j++){    
      if(j == 0 || binid[j-1] == 1){
	memcpy(Xm + mj*n, X + j*n, fbytes*n);
	//for(i = 0; i < n; i++){
	//Xm[i + mj*n] = X[i + j*n];
	//}
	mj += 1;
      }
    }

    //allocate memory for arrays required for gpuLmfit
    coeffs = (float *) calloc(pm, fbytes);
    resids = (float *) calloc(n, fbytes);
    effects = (float *) malloc(n* (*ycols) *fbytes);
    qrAux = (double *) calloc(pm, dbytes);

    if(coeffs == NULL){
      Rprintf("\nMemory for Xm coeffs to allocate\n");
      exit(1);
    }
    if(resids == NULL){
      Rprintf("\nMemory for resids failed to allocate\n");
      exit(1);
    }
    if(effects == NULL){
      Rprintf("\nMemory for effects failed to allocate\n");
      exit(1);
    }
    if(qrAux == NULL){
      Rprintf("\nMemory for qrAux failed to allocate\n");
      exit(1);
    }
    
    //required for gpuLmfit
    memcpy(effects, Y, n* (*ycols) *fbytes);

    *rank = 1; //initialize

    //required for gpuLmfit
    pivot = (int *) malloc(pm*ibytes);
    if(pivot == NULL){
      Rprintf("\nMemory for pivot failed to allocate\n");
      exit(1);
    }
    for(i = 0; i < pm; i++)
      pivot[i] = i;

    //fit the current model
    gpuLSFitF(Xm, n, pm, Y, *ycols, tol, coeffs, resids, effects, rank, pivot, qrAux);

    //compute the current model's R squared, used in marginal likelihood computation
    ssr = 0;
    for(i=0; i<n; i++)
      ssr += resids[i]*resids[i];
    sighat = ssr / (n-pm); 
    Rsq = 1 - ssr/ynorm;

    //update best set of models information if applicable
    if(*sorttype == 0){
      //save info of all models
      aics[id] = aic(n, pm, sighat);
      bics[id] = bic(n, pm, sighat);
      logmargs[id] = logmarglike(n, km, *g, Rsq);
      models[id] = id+1;
    }
    else if(*sorttype == 1){
      //save the best num_save models according to AIC
      arraysort(aics, worstid, worstscore, *num_save, 1);
      score = aic(n, pm, sighat);
      if(score < *worstscore){
	models[*worstid] = id+1;
	aics[*worstid] = score;
	bics[*worstid] = bic(n, pm, sighat);
	logmargs[*worstid] = logmarglike(n, km, *g, Rsq);
      }
    }
    else if(*sorttype == 2){
      //save the best num_save models according to BIC
      arraysort(bics, worstid, worstscore, *num_save, 1);
      score = bic(n, pm, sighat);
      if(score < *worstscore){
	models[*worstid] = id+1;
	aics[*worstid] = aic(n, pm, sighat);
	bics[*worstid] = score;
	logmargs[*worstid] = logmarglike(n, km, *g, Rsq);
      }
    }
    else{
      //save the best num_save models according to Marginal Likelihood
      arraysort(logmargs, worstid, worstscore, *num_save, 0);
      score = logmarglike(n, km, *g, Rsq);
      if(score > *worstscore){
	models[*worstid] = id+1;
	aics[*worstid] = aic(n, pm, sighat);
	bics[*worstid] = bic(n, pm, sighat);
	logmargs[*worstid] = score;
      }
    }

    //free memory use for storing information on current model
    //note: only memory whose size depends on the particular model
    free(coeffs);
    free(resids);
    free(effects);
    free(qrAux);
    free(pivot);
    free(Xm);
    
  }

  //free remaining memory
  free(rank);
  free(worstid);
  free(worstscore);

}
