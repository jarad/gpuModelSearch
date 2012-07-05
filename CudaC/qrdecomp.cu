#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<cublas.h>
#include<R.h>

#define NTHREADS 512

__device__ int getMaxIdx(int p, float * ColNorms, int id){
  int maxIdx = 0, i, bit;
  float max = ColNorms[maxIdx], tmp;

  for(i = 1; i < p; i++){
    bit =  ( id & (1 << i) ) >> 1
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

__global__ void GetColNorms(int rows, int cols, float * dX, float * ColNorms){

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

__global__ void GetQRDecompBlocked(int n, int p, float * X, float * ColNorms){
  
  float tol = 0.0001;  //float tolerance
  int blockSize = 1 << 7;  //2^7 = 128, optimal for large n according to gputools

  int blockCount;

  unsigned int id = threadIDx.x + blockIDx.x * blockDim.x;
  int varcount = 0, maxIdx, pm = 1; 
  unsigned int * binid;
  float * Xm, * pXcol = X, * pXmcol = Xm, maxNorm;
  size_t fbytes = sizeof(float);
  int rk = 0; // Local value of rank;
  int maxRank;
  float minElt; // Lowest acceptable norm under given tolerance.

  if(id > M)
    return;

  binid = (unsigned int *)calloc(p-1, sizeof(unsigned int));

  Xm = (float *) malloc(n*pm*fbytes);  //allocate memory for model matrix

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
  
  maxIdx = getMaxIdx(p, ColNorms, id);
  maxNorm = ColNorms[maxIdx];

  blockCount = (pm + blockSize - 1) / blockSize;  //set number of blocks in Xm

  if (maxNorm < tol)
    maxRank = 0; // Short-circuits the main pivot loop
  else {
    minElt = (1.0 + maxNorm) * tol;
    maxRank = rows > cols ? cols : rows; //min(rows,cols)
  }


  /*
    next steps of the algorithm (the pivoting part):

    initialize pointer to block of dQR - i.e. R. 
    note: block size must be set above first

    for each block: (  blockCount = (ncols + blockSize - 1) / blockSize;   )

    initialize int blockEnd = 0;
    
    for each unprocessed column:

  */

  



}

void gpuQRDecomp(int rows, int cols, float * X){

  
  int n = rows, p = cols, k = p - 1;
  size_t fbytes = sizeof(float);
  float * dX, * dColNorms;
  int nblocks, nthreads = NTHREADS;

  cudaMalloc( (void **)&dX, n*p*fbytes );
  cudaMalloc( (void **)&dColNorms, p*fbytes);

  cudaMemcpy( dX, X, n*p*fbytes, cudaMemcpyHostToDevice);

  nblocks = cols / nthreads;
  if(nblocks * nthreads < cols)
      nblocks++;

  getColNorms<<<nblocks, nthreads>>>(n, p, dX, dColNorms);

  
  



} 
