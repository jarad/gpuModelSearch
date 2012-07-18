#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <cublas.h>
#include <R.h>
#include "cuseful.h"


//function prototype for the function to be called from R
extern "C" void CSlmsearch(float *X, int *rows, int *cols, float *Y, int *ycols, 
			   int *g, float *aics, float *bics, float *logmargs, 
			   float *prob, float *otherprob, int *models, int *binids, 
			   int *num_save, int *sorttype );   

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

//finds the largest or smallest element in src and outputs the index and the element
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


// Rounds "length" up to the next multiple of the block length.
// from gputools
int alignBlock(int length, unsigned blockExp) {
  int blockSize = 1 << blockExp;
  return (length + blockSize - 1) & (((unsigned) -1)  << blockExp);
}


//finds column norms of a matrix, from gputools
__global__ void getColNorms(int rows, int cols, float * da, int lda, 
	float * colNorms){

	int colIndex = threadIdx.x + blockIdx.x * blockDim.x;
	float 
		sum = 0.f, term,
		* col;

	if(colIndex >= cols)
		return;

	col = da + colIndex * lda;

	// debug printing
	// printf("printing column %d\n", colIndex);
	// for(int i = 0; i < rows; i++)
	// printf("%f, ", col[i]);
	// puts("");
	// end debug printing

	for(int i = 0; i < rows; i++) {
		term = col[i];
		term *= term;
		sum += term;
	}

	// debug printing
	// printf("norm %f\n", norm);
	// end debug printing

	colNorms[colIndex] = sum;
}



// use householder xfrms and column pivoting to get the R factor of the
// QR decomp of matrix da:  Q*A*P=R, equiv A*P = Q^t * R
// (slight modification of gputools code)
__host__ void getQRDecompBlockedMod(int rows, int cols, double tol, float * dQR,
				    int blockSize, int stride, int * pivot, 
				    double * qrAux, int * rank, float *dV, float *dW, 
				    float *dT, float *du, float *dWtR, 
				    float *dColNorms){
  // Copyright 2009, Mark Seligman at Rapid Biologics, LLC.  All rights
// reserved.
//

  int
		fbytes = sizeof(float),
		rowsk = stride, // # unprocessed rows = stride - k
		colsk = cols,   // # unprocessed columns = cols - k
		maxCol = cols - 1,  // Highest candidate column:  not fixed.
		k = 0; // Number of columns processed.
  const int maxRow = rows - 1;


  checkCublasError("getQRDecompBlocked:");


    // Presets the matrix of Householder vectors, dV, to zero.
    // Padding with zeroes appears to offer better performance than
    // direct operations on submatrices:  aligned access helps
    // ensure coalescing.
    //


    checkCublasError("getQRDecompBlocked allocation:");

    // Obtains the highest valued norm in order to approximate a condition-
    // based lower bound, "minElt".
    //
    int maxIdx = cublasIsamax(cols, dColNorms, 1)-1;
    float maxNorm = cublasSnrm2(rows, dQR + stride * maxIdx, 1);
    int rk = 0; // Local value of rank;
    int maxRank;
    double minElt; // Lowest acceptable norm under given tolerance.

    if (maxNorm < tol)
      maxRank = 0; // Short-circuits the main loop
    else {
      minElt = (1.0 + maxNorm) * tol;
      maxRank = rows > cols ? cols : rows;
    }

    float * pdQRBlock = dQR;

    int blockCount = (cols + blockSize - 1) / blockSize;
    for (int bc = 0; bc < blockCount; bc++) {
      // Determines "blockEnd", which counts the number of columns remaining
      // in the upcoming block.  Swaps trivial columns with the rightmost
      // unvisited column until either a nontrivial column is found or all
      // columns have been visited.  Note that 'blockEnd <= blockSize', with
      // inequality possible only in the rightmost processed block.
      //
      // This pivoting scheme does not attempt to order columns by norm, nor
      // does it recompute norms altered by the rank-one update within the
      // upcoming block.  A higher-fidelity scheme is implemented in the non-
      // blocked form of this function.  Sapienti sat.
      //
      int blockEnd = 0;
      for (int i = k; i < k + blockSize && i < maxRank && i <= maxCol; i++) {
        float colNorm = cublasSnrm2(rows, dQR + i * stride, 1);
	while ( (colNorm < minElt) && (maxCol > i)) {
	    cublasSswap(rows, dQR + i * stride, 1, dQR + maxCol*stride, 1);
	    int tempIdx = pivot[maxCol];
	    pivot[maxCol] = pivot[i];
	    pivot[i] = tempIdx;
	    maxCol--;
	    colNorm = cublasSnrm2(rows, dQR + i * stride, 1);
        }
	if (colNorm >= minElt)
	  blockEnd++;
      }
      rk += blockEnd;
      float scales[blockSize];
      double Beta[blockSize];

      cudaMemset2D(dV, blockSize * fbytes, 0.f, blockSize * fbytes, rowsk);

      float *pdVcol = dV;
      float *pdVdiag = dV;
      float *pdQRdiag = pdQRBlock;

      for (int colIdx = 0; colIdx < blockEnd; colIdx++, pdVcol += rowsk, pdVdiag += (rowsk + 1),
           pdQRdiag += (stride + 1), k++) {

	cublasScopy(rowsk - colIdx, pdQRdiag, 1, pdVdiag, 1);

	// Builds Householder vector from maximal column just copied.
	// For now, uses slow memory transfers to modify leading element:
	//      V_1 += sign(V_1) * normV.
	//
	float v1;	  
	cublasGetVector(1, fbytes, pdVdiag, rowsk + 1, &v1, 1);
	double v1Abs = fabs(v1);
	if (k == maxRow) // The bottom row is not scaled.
	  qrAux[k] = v1Abs;
	else { // zero-valued "normV" should already have been ruled out.
	  float normV = cublasSnrm2(rowsk - colIdx, pdQRdiag, 1);
	  double recipNormV = 1.0 / normV;
	  qrAux[k] = 1.0 + v1Abs * recipNormV;
  	  scales[colIdx] = (v1 >= 0.f ? 1.f : -1.f) * recipNormV;

          // Scales leading nonzero element of vector.
	  //
	  double fac = 1.0 + normV / v1Abs;
	  cublasSscal(1, (float) fac, pdVdiag, 1);
		  
	  // Beta = -2 v^t v :  updates squared norm on host side.
	  //
	  Beta[colIdx] = -2.0 / (normV*normV + v1Abs * v1Abs * (-1.0 + fac * fac));

	  // Rank-one update of the remainder of the block, "B":
	  // u = Beta B^t v
	  //
	  cublasSgemv('T', rowsk, min(blockSize,colsk), (float) Beta[colIdx], pdQRBlock, stride, pdVcol, 1, 0.f, du, 1);
				       
          // B = B + v u^t
	  //
	  cublasSger(rowsk, min(blockSize,colsk), 1.0f, pdVcol, 1, du, 1, pdQRBlock,
	  	  stride);
	}
      }

      // If more unseen columns remain, updates the remainder of QR lying to
      // the right of the block just updated.  This must be done unless we
      // happen to have exited the inner loop without having applied any
      // Householder transformations (i.e., blockEnd == 0).
      //
      if (bc < blockCount - 1 && blockEnd > 0) {
	 // w_m = Beta (I + W V^t) v_m, where the unsubscripted matrices
	 // refer to those built at step 'i-1', having 'i' columns.
	 //
	 // w_i = Beta v_i
	 //
	 //  T = V^t V
	 //
         cublasSsyrk('U', 'T', blockSize, rowsk, 1.f, dV, rowsk, 0.f, dT, blockSize);

	 float *pdTcol = dT;
	 float *pdWcol = dW;
	 pdVcol = dV;
	 for (int m = 0; m < blockSize; m++, pdWcol += rowsk, pdVcol += rowsk, pdTcol += blockSize) {
	   cublasScopy(rowsk, pdVcol, 1, pdWcol, 1);
           cublasSscal(rowsk, Beta[m], pdWcol, 1);
	   // w_m = w_m + Beta W T(.,m)
	   //
	   if (m > 0) {
	     cublasSgemv('N', rowsk, m, Beta[m], dW, rowsk, pdTcol, 1, 1.f, pdWcol, 1);
	   }
	 }

	 // Updates R, beginning at current diagonal by:
	 //   R = (I_m + V W^t) R = R + V (W^t R)
	 //

	 // WtR = W^t R
	 //
	 cublasSgemm('T','N', blockSize, colsk - blockSize, rowsk, 1.f, dW, rowsk, 
	    pdQRBlock + blockSize * stride, stride, 0.f, dWtR, blockSize);

	 // R = V WtR + R
	 //
	 cublasSgemm('N', 'N', rowsk, colsk - blockSize, blockSize, 1.f, dV,
	       rowsk, dWtR, blockSize, 1.f, pdQRBlock+ blockSize * stride, stride);
     }

      // Flushes scaled Householder vectors to the subdiagonals of dQR,
      // 'blockSize'-many at a time.  The only time a smaller number are
      // sent occurs when a partial block remains at the right end.
      //
      pdVdiag = dV;
      pdQRdiag = pdQRBlock;
      for (int l = 0; l < blockEnd; l++, pdVdiag += (rowsk + 1), pdQRdiag += (stride + 1)) {
	cublasSscal(rowsk - (l + 1), scales[l], pdVdiag + 1, 1);
	cublasScopy(rowsk - (l + 1), pdVdiag + 1, 1, pdQRdiag + 1, 1);
      }

      pdQRBlock += blockSize * (stride + 1);
      colsk -= blockSize;
      rowsk -= blockSize;
    }

    *rank = rk;	
    checkCublasError("getQRDecompBlocked, postblock:");
    checkCudaError("getQRDecompBlocked:");
     
    // dQR now contains the upper-triangular portion of the factorization,
    // R.
    // dV is lower-triangular, and contains the Householder vectors, from
    // which the Q portion can be derived.  An adjusted form of the
    // diagonal is saved in qrAux, while the sub-diagonal portion is
    // written onto QR.


    checkCublasError("getQRDecompBlocked, freed memory:");
}


// Fills in the coefficients, residuals and effects matrices.
// A slight tweak of a gputools function
__host__ void getCREmod(float *dQR, int rows, int cols, int stride, int rank, 
			double *qrAux, int yCols, float *coeffs, float *resids, 
			float *effects, float *dDiags, float *dCoeffs, float *dResids,
			float *dEffects){
	const int
		fbytes = sizeof(float);
        // Used by effects, residual computations.
	//
	int maxIdx = min(rank, rows - 1);

	float
	  * diags = Calloc(rank * fbytes, float);


	// Temporarily swaps diagonals with qrAux.

	cublasScopy(rank, dQR, stride + 1, dDiags, 1);
	cublasGetVector(rank, fbytes, dDiags, 1, diags, 1);

	float *qrAuxFloat = Calloc(maxIdx * fbytes, float);
	for (int i = 0; i < maxIdx; i++)
	  qrAuxFloat[i] = qrAux[i];
	cublasSetVector(maxIdx, fbytes, qrAuxFloat, 1, dQR, stride + 1);
	Free(qrAuxFloat);


	// Computes the effects matrix, intialized by caller to Y.

	float
		* pEffects = dEffects;

	for (int i = 0; i < yCols; i++, pEffects += rows) {
		float
			* pQR = dQR;

		for (int k = 0; k < maxIdx; k++, pQR += (stride + 1)) {
			double
				t = cublasSdot(rows - k, pQR, 1, pEffects +  k, 1);

			t *= -1.0 / qrAux[k];
			cublasSaxpy(rows - k, t, pQR, 1, pEffects + k, 1);
		}
	}

	
	// Computes the residuals matrix, initialized by caller to zero.
	// If not of full row rank, presets the remaining rows to those from
	// effects.
	if(rank < rows) {
		for(int i = 0; i < yCols; i++) {
			cublasScopy(rows - rank,  dEffects + i*rows + rank, 1,
				dResids + i*rows + rank, 1);
		}
	}

	float
		* pResids = dResids;


	for (int i = 0; i < yCols; i++, pResids += rows) {
		for (int k = maxIdx - 1; k >= 0; k--) {
			double
				t = -(1.0 / qrAux[k])
					* cublasSdot(rows - k, dQR + k*stride + k, 1, pResids + k, 1);

			cublasSaxpy(rows -k, t, dQR + k*stride + k, 1, pResids + k, 1);
		}
	}
	cublasScopy(maxIdx, dDiags, 1, dQR, stride + 1);

	// Computes the coefficients matrix, initialized by caller to zero.

	float
		* pCoeffs = dCoeffs;

	for(int i = 0; i < yCols; i++, pCoeffs += cols) {
		cublasScopy(rank, dEffects + i*rows, 1, pCoeffs, 1);

		float t;
		for(int k = rank - 1; k > 0; k--) {
		        cublasSscal(1, 1.f / diags[k], pCoeffs + k, 1);
			cublasGetVector(1, fbytes, pCoeffs + k, 1, &t, 1);
			cublasSaxpy(k, -t, dQR + k*stride, 1, pCoeffs, 1);
		}
		cublasSscal(1, 1.f / diags[0], pCoeffs, 1);
	}
	Free(diags);

	cublasGetMatrix(cols, yCols, fbytes, dCoeffs, cols, coeffs, cols);
	cublasGetMatrix(rows, yCols, fbytes, dResids, rows, resids, rows);
	cublasGetMatrix(rows, yCols, fbytes, dEffects, rows, effects, rows);

	


}



//performs model search
void CSlmsearch(float *X, int *rows, int *cols, float *Y, int *ycols, int *g, 
		float *aics, float *bics, float *logmargs, float *prob, 
		float *otherprob, int *models, int *binids, int *num_save, 
		int *sorttype ){

  int p = *cols, k = p - 1, M = 1 << k, n = *rows;
  int id, km, pm, mj;
  int i, j;
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
  float totalprob;
  float *dX, *dXm, *dY;
  const unsigned blockExp = 7; // Gives blockSize = 2^7 = 128.
  int stride = alignBlock(n, blockExp);
  int blockSize = 1 << blockExp;
  float *dV, *dW, *dT, *du, *dWtR, *dColNorms, *dFColNorms;
  int nthreads = 512, nblocks;
  float *dResids, *dEffects, *dCoeffs, *dDiags;

  //sets the number of blocks to use in column norm algorithm
  nblocks = p / nthreads;
  if(nblocks * nthreads < p)
    nblocks++;

  cublasInit();

  cublasAlloc(stride * p, fbytes, (void **)&dX);
  cublasAlloc(stride * p, fbytes, (void **)&dXm);
  cublasAlloc(stride * blockSize, fbytes, (void**) &dV);
  cublasAlloc(stride * blockSize, fbytes, (void**) &dW);
  cublasAlloc(blockSize * blockSize, fbytes, (void**) &dT);
  cublasAlloc(stride, fbytes, (void**) &du);
  cublasAlloc(blockSize * (p - blockSize), fbytes, (void**) &dWtR);
  cublasAlloc(p, fbytes, (void**) &dColNorms);
  cublasAlloc(p, fbytes, (void**) &dFColNorms);
  cublasAlloc(n * (*ycols), fbytes, (void**) &dResids);
  cublasAlloc(n * (*ycols), fbytes, (void**) &dEffects);
  cublasAlloc(p * (*ycols), fbytes, (void **) &dCoeffs);
  cublasAlloc(p, fbytes, (void**) &dDiags);
  cublasAlloc(n * (*ycols), fbytes, (void**) &dY);
  
  // This is overkill:  just need to zero the padding.
  // (copy full matrix X to device)
  cudaMemset2D(dX, p * fbytes, 0.f, p * fbytes, stride); 
  cublasSetMatrix(n, p, fbytes, X, n, dX, stride);

  //copy Y to device
  cublasSetMatrix(n, (*ycols), fbytes, Y, n, dY, n);

  //calculate column norms of full matrix
  getColNorms<<<nblocks, nthreads>>>(n, p, dX, stride, dFColNorms);
    
  rank = (int *) malloc(ibytes);
  worstid = (int *) malloc(ibytes);
  worstscore = (float *) malloc(fbytes);

  //allocate memory for arrays required for getQR
  coeffs  = (float *)  malloc(p * (*ycols) * fbytes);
  resids  = (float *)  malloc(n * fbytes);
  effects = (float *)  malloc(n * (*ycols) *  fbytes);
  qrAux   = (double *) malloc(p * dbytes);
  pivot   = (int *)    malloc(p * ibytes);

  if(pivot == NULL){
    Rprintf("\nMemory for pivot failed to allocate\n");
    exit(1);
  }
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

    //copy the relevant columns of the full model matrix into the current model's
    //model matrix (Xm). Note: R stores matrices/arrays in column major format.
    //Also calculates pm=ncol(Xm) and extracts appropriate column norms.
    mj = 0; //keeps track of column numbers of current model matrix
    pm = 0; //keeps track of number of variables in model matrix

    for(j = 0; j < p; j++){  
      if(j == 0){
	//copy intercept column
	cublasScopy(stride, dX + j*stride, 1, dXm + mj*stride, 1);
	cublasScopy(1, dFColNorms + j, 1, dColNorms + mj, 1);
	pm += 1;
	mj += 1;
      }//end if
      else{
	//copy covariate columns
	bit = ( id & ( 1 << (j - 1) ) )  >> (j - 1) ;
	binid[k - j] = bit;
	if(bit == 1){
	  cublasScopy(stride, dX + j*stride, 1, dXm + mj*stride, 1);
	  cublasScopy(1, dFColNorms + j, 1, dColNorms + mj, 1);
	  pm += 1;
	  mj += 1;
	}//end if

      }//end else

    }//end for loop
    
    km = pm - 1;
    
    //initialize
    *rank = 1;
    for(i = 0; i < pm; i++)
      pivot[i] = i;

    //initialize device coefficients and resids to 0, effects to Y
    cudaMemset(dCoeffs, 0.f, pm * (*ycols) * fbytes);
    cudaMemset(dResids, 0.f, n  * (*ycols) * fbytes);
    cublasScopy(n*(*ycols), dY, 1, dEffects, 1);

    //resetting to 0
    for(i = 0; i < pm; i++)
      qrAux[i] = 0;

    // On return we have dXm in pivoted, packed QR form.
    //
    getQRDecompBlockedMod(n, pm, tol, dXm, blockSize, stride, pivot, qrAux, rank, 
			  dV, dW, dT, du, dWtR, dColNorms);
    
    //get coefficients, residuals, and effects for current model using QR decomp
    if(*rank > 0)
      getCREmod(dXm, n, pm, stride, *rank, qrAux, *ycols, coeffs, resids, effects,
	dDiags, dCoeffs, dResids, dEffects);
    else // Residuals copied from Y.
      memcpy(resids, Y, n * (*ycols) * fbytes);


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
      if(id == 1){
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
  free(rank);
  free(worstid);
  free(worstscore);
  cublasFree(dX);
  cublasFree(dXm);
  cublasFree(dW);
  cublasFree(dV);
  cublasFree(dT);
  cublasFree(du);
  cublasFree(dWtR);
  cublasFree(dColNorms);
  cublasFree(dFColNorms);
  cublasFree(dResids);
  cublasFree(dEffects);
  cublasFree(dCoeffs); 
  cublasFree(dY);
  free(coeffs);
  free(resids);
  free(effects);
  free(qrAux);
  free(pivot);
  
  
  cublasShutdown();

}

