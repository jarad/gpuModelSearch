#include<stdio.h>
#include<stdlib.h>
#include<string.h>
//#include<cublas.h>
//#include<R.h>

#define NTHREADS 512


//computes u = constant * t(X) %*% v
__device__ void cXtv(float con, int rows, int cols, float * X, int ldX, float * v, 
		     float * u){
  int i,k;
  float sum;

  for(k = 0; k < cols; k++){
    sum = 0;
    for(i = 0; i < rows; i++){
      sum += X[i + k*ldX] * v[i];
    }
    u[k] = con*sum;
  }
}

//computes X = X + v %*% u'
__device__ void Xpvut(int rows, int cols, float * X, int ldX, float *v, float *u){
  int i, j;
  
  for(i=0; i<rows; i++)
    for(j=0; j<cols; j++)
      X[i + j*ldX] += v[i] * u[j];
}


//calculates T = V'V where V is upper triangular
__device__ void uVtV(int rows, int cols, float * V, float * T){
  int i, j, k;
  float sum;
  for(i=0; i<cols; i++){
    for(j=0; j<cols; j++){
      sum = 0;
      for(k=0; k<min(i,j); k++)
	sum += V[k + i*rows] * V[k + j*rows];
      T[i + j*cols] = sum;
    }
  }
}


//Rounds "length" up the the next multiple of block length
__device__ int alignBlock(int length, unsigned blockExp){
  int blockSize = 1 << blockExp;
  int out = (length + blockSize - 1) & ( ( (unsigned) - 1) << blockExp);
  return out;
}

__device__ void colswap(float *X, int n, int p, int i, int j){
  float tmp;
  int k;
  float * Xi = X + i*n;
  float * Xj = X + i*j;
  if(i <= p && j <= p){
    for(k = 0; k < n; k++){
      tmp = Xi[k];
      Xi[k] = Xj[k];
      Xj[k] = tmp;
    }
  }
}

__device__ int getMaxIdx(int p, float * ColNorms, int id){
  int maxIdx = 0, i, bit;
  float max = ColNorms[maxIdx], tmp;

  for(i = 1; i < p; i++){
    bit =  ( id & (1 << i) ) >> 1;
    if( bit == 1){
      tmp = ColNorms[i];
      if(tmp > max){
	max = tmp;
	maxIdx = i;
      }
    }
  }

  return maxIdx;
}

__device__ int idx(int id, int p, int i){
  int j, colcount=0, out;

  for(j=0; j<p; j++){

    if(colcount == i)
      out = j;
    
    if( ( id & ( 1 << j ) ) >> j == 1)
      colcount += 1;
  }

  return out;
}

__global__ void getColNorms(int rows, int cols, float * dX, float * ColNorms){

	int colIndex = threadIdx.x + blockIdx.x * blockDim.x;
	float sum = 0.f, term, * col;

	if(colIndex >= cols)
	  return;

	col = dX + colIndex * rows;

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

	ColNorms[colIndex] = sqrtf(sum);
}

__global__ void GetQRDecompBlocked(int n, int p, float * X, float * ColNorms, 
				   unsigned int blockSize, int stride, float * scales, 
				   float * Beta, float * qrAux, int * pivot, 
				   float * WtR, int * rank, float * V, float * W,
				   float * T, float * u){
  
  float tol = 0.0001;  //float tolerance
  int blockCount, blockEnd = 0, bc, i, kp;  //kp is # cols processed

  //unsigned int id = threadIDx.x + blockIDx.x * blockDim.x;
  int varcount = 0, maxIdx, pm = 1; 
  //unsigned int * binid;
  float * Xm; //model id's model matrix; will eventually contain R in QR decomp
  //float * pXcol = X; //pointer to column of X
  float * pXmcol; //pointer to column of Xm
  float maxNorm;
  float * pXmBlock; //pointer to block of Xm 
  size_t fbytes = sizeof(float);
  int rk = 0; // Local value of rank;
  int maxRank;  //maximum possible rank of the model matrix
  float minElt; // Lowest acceptable norm under given tolerance.
  int maxCol; //this will be initialized to pm - 1, highest candidate column, not fixed
  int maxRow = n - 1;
  int rowskp = stride; // # unprocessed rows = stride - kp
  int colskp;   // # unprocessed cols = cols - kp
  float *pVcol;             //pointer to column of V
  float *pVdiag;            //pointer to diagonal element of V
  float *pXmdiag;    //pointer to diagonal element of dQRBlock
  int colIdx;
  int a,b;
  int l;


  //if(id > M)
  //  return;

  //binid = (unsigned int *)calloc(p-1, sizeof(unsigned int));

  /*
  //copies the relevant columns of X into Xm
  do{
    int jump = 0;
    unsigned int bit = 0;

    memcpy(pXmcol, pXcol, n*fbytes);  //copy the current column

    //any time we copy any column but the first, set binid[col] = 1
    if(varcount > 0){ 
      binid[p - varcount] = 1; //actual index is (p - 1) - (varcount - 1)
      pm += 1;
    }

    pXmcol += n;  //update pointer on Xm to the next column

    //bitwise operations on id to find proper columns
    //look for the next used column of X
    do{
      jump += 1;  //increment number of columns of X to move forward
      varcount += 1; //increment number of total vars of X checked
      //note: each column of X but the first represents a variable

      //set bit as the k'th bit of id (1st bit is far right)
      bit = ( id & ( 1 << varcount ) ) >> varcount; 

      //alternative code: declare tmpid above first.
      ////bit = tmpid & 1;  //bit = 1 if last bit of tmpid is 1, 0 otherwise
      ////tmpid = tmpid >> 1;  //bitshift tmpid right - lopping off last bit
    }
    while(bit == 0 && varcount < p);

       
    pXcol += n*jump; //update pointer on X to the next used column
  }
  while(varcount < p );  
  //if varcount = k ( = p - 1 ), we still have to copy this last column, so the
  //loop must iterate again
  */

  pm = p;
  
  //initialize variables that depend on pm and Xm
  //Xm = (float *) calloc(stride*pm,fbytes);  //allocate memory for model matrix
  Xm = X;
  pXmcol = Xm;
  pXmBlock = Xm; 

  maxCol = pm - 1; //maximum possible column

  colskp = pm;
  maxIdx = getMaxIdx(p, ColNorms, (1 << pm) - 1);
		     //id);
  maxNorm = ColNorms[maxIdx];

  if (maxNorm < tol){
    maxRank = 0; // Short-circuits the main pivot loop
  }
  else{
    minElt = (1.0 + maxNorm) * tol;
    maxRank = n > pm ? pm : n; //min(n,pm)
  }

  blockCount = (pm + blockSize - 1) / blockSize;
  
  for(bc = 0; bc < blockCount; bc++){
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

    blockEnd = 0;
    
    for(i = kp; i < kp + blockSize && i < maxRank && i <= maxCol; i++){
      float colNorm = ColNorms[i];
	//ColNorms[idx(id, p, i)];

      while( (colNorm < minElt) && (maxCol > i) ){
	//while colNorm is smaller than the lowest acceptable norm, given tolerance
	//and while the maximum column possible is larger than the current column

	//keep track of column swap in pivot vector
	int tempIdx = pivot[maxCol];
	
	//exchange columns of Xm: current column and maxCol
	colswap(Xm, stride, pm, i, maxCol);
	
	pivot[maxCol] = pivot[i];
	pivot[i] = tempIdx;
	    
	maxCol--;

	//set new colNorm as norm of max column, now in i'th spot
	colNorm = ColNorms[i];
	  //ColNorms[idx(id, p, i)];
      }

      if (colNorm >= minElt)
	blockEnd++;
       
    }

    rk += blockEnd;

    //set V to 0 everywhere
    for(a=0; a<blockSize; a++)
      for(b=0; b<rowskp; b++)
	*(V + b + a*rowskp) = 0;

    
    pVcol = V;             //pointer to column of V
    pVdiag = V;            //pointer to diagonal element of V
    pXmdiag = pXmBlock;    //pointer to diagonal element of XmBlock
    
    for(colIdx = 0; colIdx < blockEnd; colIdx++){
      float v1;
      float v1Abs;

      memcpy(pVdiag, pXmdiag, (rowskp - colIdx)*fbytes);
      
      // Builds Householder vector from maximal column just copied.
      // For now, uses slow memory transfers to modify leading element:
      //      V_1 += sign(V_1) * normV.

      v1 = *pVdiag;
      v1Abs = fabs(v1);

      if(kp == maxRow)//botom row is not scaled
	qrAux[kp] = v1Abs;
      else{ // zero-valued "normV" should already have been ruled out.
	int d;
	float normV = 0;
	float recipNormV;
	float fac;
		
	for(d=0; d < rowskp - colIdx; d++){
	  float diag = pVdiag[d];
	  normV += diag*diag;
	}

	normV = sqrtf(normV);
	recipNormV = 1.0 / normV;
	qrAux[kp] = 1.0 + v1Abs * recipNormV;
	scales[colIdx] = (v1 >= 0.f ? 1.f : -1.f) * recipNormV;
	
	// Scales leading nonzero element of vector.
	fac = 1.0 + normV / v1Abs;
	*pVdiag *= fac;

	Beta[colIdx] = -2.0 / (normV * normV + v1Abs * (-1.0 + fac * fac));

	// u = Beta[colIdx] * t(pXmBlock) %*% pVcol
	// rows of pXmBlock = rowskp
	// cols of pXmBlock = min(blockSize, colskp)
	// leading dimension of pXmBlock = stride
	// elements of pVcol = rowskp
	// elements of u = min(blockSize, colskp)
	cXtv(Beta[colIdx], rowskp, min(blockSize, colskp), pXmBlock, stride, pVcol, u);
	
	//pXmBlock = pXmBlock + pVcol %*% u'
	//rows of pXmBlock = rowskp
	//cols of pXmBlock = min(blockSize, colsk)
	//elements of pVcol = rowskp
	//elements of u = min(blockSize, colsk)
	//leading dim of pXmBlock = stride
	Xpvut(rowskp, min(blockSize, colskp), pXmBlock, stride, pVcol, u);
	
      }

      pVcol += rowskp; 
      pVdiag += (rowskp + 1); 
      pXmdiag += (stride + 1); 
      kp++;
    }

    
    // If more unseen columns remain, updates the remainder of QR lying to
    // the right of the block just updated.  This must be done unless we
    // happen to have exited the inner loop without having applied any
    // Householder transformations (i.e., blockEnd == 0).
    if (bc < blockCount - 1 && blockEnd > 0) {
      // w_m = Beta (I + W V^t) v_m, where the unsubscripted matrices
      // refer to those built at step 'i-1', having 'i' columns.
      //
      // w_i = Beta v_i
      //

      
      float *pTcol = T;  //pointer to columns of dT
      float *pWcol = W;  //pointer to columns of dW
      int m, a, b, c;
      pVcol = V;  //pointer to columns of dV
      
      //  T = V^t V
      //// T = dT
      //// V = dV
      //rows of V' = blockSize
      //cols of V' = rowsk
      //leading dim of V = rowsk
      //leading dim of T = rowsk
      //rows/cols of T = blockSize
      uVtV(rowskp, blockSize, V, T);

      for (m = 0; m < blockSize; m++, pWcol += rowskp, pVcol += rowskp, 
	     pTcol += blockSize) {
	//for m=0,1,..., blockSize - 1
	//pdWcol += rowsk
	//pdVcol += rowsk
	//pdTcol += blockSize
	int a;

	for(a = 0; a<rowskp; a++)
	  pWcol[a] = Beta[m]*pVcol[a];
	
	// w_m = w_m + Beta W T(.,m)
	//// w_m = pWcol
	//// Beta = Beta[m]
	//// T(.,m) = pTcol
	//// W = dW
	if (m > 0) {
	  int a, b;
	  for(a = 0; a < rowskp; a++){
	    float sum = 0;
	    for(b = 0; b < m; b++){
	      sum += W[b + a*m] * pTcol[b];
	    }
	    pWcol[a] += sum;
	  }
	}
      }

      
      // Updates R, beginning at current diagonal by:
      //   R = (I_m + V W^t) R = R + V (W^t R)
      //
      
      // WtR = W^t R
      //// WtR = dWtR
      //// W = dW
      //W is rowskp by blockSize
      //pXmBlock + blockSize*stride is rowskp by colskp - blockSize
      //WtR is blockSize by colsk - blockSize
      for(a = 0; a < blockSize; a++){
	for(b = 0; b < colskp - blockSize; b++){
	  float sum = 0;
	  for(c = 0; c < rowskp; c++)
	    sum += W[c + rowskp*a] * *(pXmBlock + blockSize*stride + c + b*rowskp);
	  WtR[a + b*blockSize] = sum;
	}
      }

      // R = V WtR + R
      //// V = dV 
      //// WtR = dWtR
      //// R = pdQRBlock + blockSize*stride
      // V is rowskp by blockSize
      //WtR is blockSize by colskp - blockSize
      //pXmBlock + blockSize*stride is rowskp by colskp - blockSize
      for(a = 0; a < rowskp; a++){
	for(b = 0; b < colskp - blockSize; b++){
	  float sum = 0;
	  for(c = 0; c < blockSize; c++){
	    sum += V[a + c*rowskp] * WtR[c + b*blockSize];
	  }
	  *(pXmBlock + blockSize*stride + a + b*rowskp) += sum;
	}
      }

    }
      
    // Flushes scaled Householder vectors to the subdiagonals of dQR,
    // 'blockSize'-many at a time.  The only time a smaller number are
    // sent occurs when a partial block remains at the right end.
    //
    pVdiag = V;  //pointer to diagonal of V
    pXmdiag = pXmBlock;  //pointer to diagonal of pXmBlock / XmBlock

    
    for (l = 0; l < blockEnd; l++, pVdiag += (rowskp + 1),
	   pXmdiag += (stride + 1)) {
      //for l = 0, 1, ..., blockEnd - 1
      //// pVdiag += (rowskp + 1)
      //// pXmdiag += (stride +1)
      
      int a;

      //pVdiag + 1 = scales[l]* (pVdiag + 1)
      //pXmdiag + 1 = pVdiag + 1
      for(a = 0; a < rowskp - l - 1; a++){
	*(pVdiag + 1 + a) *= scales[l];
	*(pXmdiag + 1 + a) = *(pVdiag + 1 + a);
      }

    }
    
    pXmBlock += blockSize * (stride + 1);
    colskp -= blockSize;
    rowskp -= blockSize;
    
    //end main loop

  }

  
  //set rank
  *rank = rk;	
       
  // Xm now contains the upper-triangular portion of the factorization,
  // R.
  // V is lower-triangular, and contains the Householder vectors, from
  // which the Q portion can be derived.  An adjusted form of the
  // diagonal is saved in qrAux, while the sub-diagonal portion is
  // written onto QR.


}

int main(void){
  int n = 10, p = 5;
  //, k = p - 1;
  int i, j;
  size_t fbytes = sizeof(float);
  float * X = (float *)malloc(n*p*fbytes);
  float * dX, * dColNorms;
  int nblocks, nthreads = NTHREADS;
  const unsigned blockExp = 7;
  unsigned int blockSize = 1 << blockExp;
  int stride = (n + blockSize - 1) & (((unsigned) -1) << blockExp);

  float * dscales;
  float * dBeta;
  float * dqrAux;
  int * dpivot;
  float * dWtR;
  int * drank;
  float * dV;
  float * dW; 
  float * dT; 
  float * du; 

  for(i=0; i<n; i++)
    for(j=0; j<p; j++)
      X[i + j*n] = i + j;

  printf("\n X \n");

  for(i=0; i<n; i++){
    for(j=0; j<p; j++){
      printf("%.f ", X[i + j*n]);
    }
    printf("\n");
  }
    
  
  cudaMalloc( (void **)&dX, n*p*fbytes );
  cudaMalloc( (void **)&dColNorms, p*fbytes);
  cudaMalloc( (void **)&dscales, blockSize*fbytes);
  cudaMalloc( (void **)&dBeta, blockSize*fbytes);
  cudaMalloc( (void **)&dqrAux, p*fbytes);
  cudaMalloc( (void **)&dpivot, p*sizeof(int));
  cudaMalloc( (void **)&dWtR, blockSize*(p - blockSize)*fbytes);
  cudaMalloc( (void **)&drank, sizeof(int));
  cudaMalloc( (void **)&dV, blockSize*stride*fbytes);
  cudaMalloc( (void **)&dW, blockSize*stride*fbytes);
  cudaMalloc( (void **)&dT, blockSize*blockSize*fbytes);
  cudaMalloc( (void **)&du, stride*fbytes);


  cudaMemcpy( dX, X, n*p*fbytes, cudaMemcpyHostToDevice);

  nblocks = p / nthreads;
  if(nblocks * nthreads < p)
      nblocks++;

  getColNorms<<<nblocks, nthreads>>>(n, p, dX, dColNorms);

  GetQRDecompBlocked<<<1,1>>>(n, p, dX, dColNorms, blockSize, stride, dscales, dBeta,
			      dqrAux, dpivot, dWtR, drank, dV, dW, dT, du);

  cudaMemcpy( X, dX, n*p*fbytes, cudaMemcpyDeviceToHost);

  printf("\n R \n");

  for(i=0; i<n; i++){
    for(j=0; j<p; j++){
      printf("%.f ", X[i + j*n]);
    }
    printf("\n");
  }

  cudaFree(dX);
  cudaFree(dColNorms);
  free(X);
  cudaFree(dscales);
  cudaFree(dBeta);
  cudaFree(dqrAux);
  cudaFree(dpivot);
  cudaFree(dWtR);
  cudaFree(drank);
  cudaFree(dV);
  cudaFree(dW);
  cudaFree(dT);
  cudaFree(du);

  return 0;

} 
