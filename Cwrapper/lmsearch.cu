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
extern "C" void Clmsearch(float *X, int *rows, int *cols, float *Y, int *ycols, 
			  int *g, float *aics, float *bics, float *logmargs, 
			  double *prob, double *otherprob, int *models, int *binids,
			  int *num_save, int *sorttype);   

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
void Clmsearch(float *X, int *rows, int *cols, float *Y, int *ycols, int *g, 
	       float *aics, float *bics, float *logmargs, double *prob, 
	       double *otherprob, int *models, int *binids, int *num_save, 
	       int *sorttype){

  int p = *cols, k = p - 1, M = 1 << k, n = *rows;
  int id, km, pm, mj;
  int i, j;
  float *Xm; //pointer to model matrix of current model
  size_t fbytes = sizeof(float), ibytes = sizeof(int), dbytes = sizeof(double);
  float *coeffs, *resids, *effects; //pointers to arrays, needed for gpuLmfit
  int *rank, *pivot; 
  double *qrAux;
  double tol = 0.0001; //tolerance of the fit, gputools only supports float
  float score, *worstscore; 
  float sighat, ynorm=0, ssr, ymn=0, Rsq;
  int *worstid;
  int bit, binid[k];
  float a, b, lml, maxlml;
  double totalprob;

  rank       = (int *)  malloc(ibytes);
  worstid    = (int *)  malloc(ibytes);
  worstscore = (float *) malloc(fbytes);
  
  //required for gpuLmfit
  pivot = (int *) malloc(p * ibytes);
  
  //allocate memory for model matrix of current model
  //note: will be overwritten repeatedly in loop
  Xm = (float *) malloc( n * p * fbytes );
  
  //allocate memory for arrays required for gpuLmfit
  //note: will be overwritten repeatedly in loop
  coeffs  = (float *)  malloc(p * (*ycols) * fbytes);
  resids  = (float *)  malloc(n * (*ycols) * fbytes);
  effects = (float *)  malloc(n * (*ycols) * fbytes);
  qrAux   = (double *) malloc(p * dbytes);
  
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
  if(Xm == NULL){
    Rprintf("\nMemory for Xm failed to allocate\n");
    exit(1);
  }
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
  if(pivot == NULL){
    Rprintf("\nMemory for pivot failed to allocate\n");
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

    mj = 0; //keeps track of column numbers of current model matrix
    pm = 0;
    
    //copy the relevant columns of the full model matrix into the current model's
    //model matrix. Note: R stores matrices/arrays in column major format.
    for(j = 0; j < p; j++){    
      if(j == 0){
	memcpy(Xm + mj*n, X + j*n, fbytes*n);
	mj += 1;
	pm += 1;
      }//end if
      else{
	bit = (id & ( 1 << (j - 1) ) ) >> (j - 1);
	binid[k - j] = bit;
	if(bit == 1){
	  memcpy(Xm + mj*n, X + j*n, fbytes*n);
	  mj += 1;
	  pm += 1;
	}//end if
      }//end else
    }//end loop
    km = pm - 1;

    //initialization required for gpuLmfit
    //note: not using memset to avoid possible unwanted behavior
    memcpy(effects, Y, n * (*ycols) *fbytes);
    for(i = 0; i < n; i++){
      if(i < pm){
	coeffs[i] = 0;
	qrAux[i] = 0;
      }
      else if(i < pm * (*ycols) ){
	coeffs[i] = 0;
      }
				  
      resids[i] = 0;
    }

    *rank = 1; //initialize
  
    for(i = 0; i < pm; i++)
      pivot[i] = i;

    //fit the current model
    gpuLSFitF(Xm, n, pm, Y, *ycols, tol, coeffs, resids, effects, rank, pivot, qrAux);

    //compute the current model's R squared, used in marginal likelihood computation
    ssr = 0;
    for(i = 0; i < n * (*ycols); i++)
      ssr += resids[i]*resids[i];
    sighat = ssr / (n-pm); 
    Rsq = 1 - ssr/ynorm;

    a = aic(n, pm, sighat);
    b = bic(n, pm, sighat);
    lml = logmarglike(n, km, *g, Rsq);

    if(lml > maxlml){
      if(id == 0){
        maxlml = lml;
        totalprob = 1;
      }
      else{
        totalprob = totalprob * exp(maxlml - lml) + 1;
        maxlml = lml;
      }
    }
    else{
      totalprob += exp(lml - maxlml);
    }

    //update best set of models information if applicable
    if(*sorttype == 0){
      //save info of all models
      aics[id] = a;
      bics[id] = b;
      logmargs[id] = lml;
      models[id] = id+1;

      //copy binary ID of current model
      for(i = 0; i < k; i++){
	binids[id + i * (*num_save)] = binid[i];
      }


    }
    else{ 

      if(*sorttype == 1){
	//save the best num_save models according to AIC
	arraysort(aics, worstid, worstscore, *num_save, 1);
	score = a;
      }
      else if(*sorttype == 2){
	//save the best num_save models according to BIC
	arraysort(bics, worstid, worstscore, *num_save, 1);
	score = b;
      }
      else{
	//save the best num_save models according to Marginal Likelihood
	arraysort(logmargs, worstid, worstscore, *num_save, 0);
	score = lml;
      }

      if(score < *worstscore){
	models[*worstid] = id+1;
	aics[*worstid] = a;
	bics[*worstid] = b;
	logmargs[*worstid] = lml;

	//copy binary ID of current model
	for(i = 0; i < k; i++){
	  binids[*worstid + i * (*num_save)] = binid[i];
	}

      }

    }
    
  }

  //calculate probabilities
  for(i = 0; i < *num_save; i++){
    prob[i] = exp(logmargs[i] - maxlml) / totalprob;
    *otherprob += prob[i];
  }
  *otherprob = 1 - *otherprob;

  //free memory
  free(Xm);
  free(rank);
  free(worstid);
  free(worstscore);
  free(coeffs);
  free(resids);
  free(effects);
  free(qrAux);
  free(pivot);


}
