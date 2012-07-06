cXtv <- function(con, rows, cols, X, v){
  u <- rep(0,cols)
  for(k in 1:cols){
    sum <- 0
    for(i in 1:rows){
      sum <- sum + v[i]*X[i,k]
    }
    u[k] <- sum*con
  }
  return(u)
}

Xpvut <- function(rows, cols, X, v, u){
  out <- X
  for(i in 1:rows)
    for(j in 1:cols)
      out[i,j] <- v[i]*u[j] + out[i,j]

  return(out)
}

uVtV <- function(rows, cols, V){
  T <- matrix(0,cols,cols)
  for(i in 1:cols){
    for(j in 1:cols){
      sum <-  0
      for(k in 1:min(i,j)){
        sum <-  sum + V[k,i]*V[k,j]
      }
      T[i,j] <- sum
    }
  }
  return(T)
}

colswap <- function(X, n, p, i, j){
  out <- X
  for(k in 1:n){
    tmp <-  out[k,i]
    out[k,i] <- out[k,j]
    out[k,j] <- tmp
  }
  return(out)
}

getMaxIdx <- function(p, ColNorms){
  out <- 0
  for(i in 1:p){
    tmp <- ColNorms[i]
    if(tmp > out)
      out <-  i
  }
  return(out)
}

getColNorms <- function(rows, cols, dX){
  ColNorms <- rep(0,cols)
  for(i in 1:cols){
    sum <- 0
    for(j in 1:rows){
      term = dX[j,i]
      sum <-  sum + term*term
    }
    ColNorms[i] <- sqrt(sum)
  }
  return(ColNorms)
}

GetQRDecompBlocked <- function(n, p, X, ColNorms, blockSize, scales, Beta,
                               qrAux, pivot, WtR, rank, V, W, T, u){
  tol <- 0.0001
  blockEnd <- 0
  rk <- 0
  rowskp <- n
  maxCol <- p 
  colskp <- p
  maxIdx <- getMaxIdx(p, ColNorms)
  maxNorm <- ColNorms[maxIdx]
  maxRow <- n
  kp <- 0

  pXmcol <- c(1,1)
  pXmBlock <- c(1,1)
  

  if(maxNorm < tol)
    maxRank <- 0
  else{
    minElt <- tol*(1.0+maxNorm)
    maxRank <- min(n,p)
  }

  blockCount <- (p + blockSize - 1) %/% blockSize

  for(bc in 1:blockCount){
    blockEnd <- 0

    for(i in (kp + 1):min(kp+blockSize, maxRank, maxCol+1) ){
      colNorm <- ColNorms[i]
      while( (colNorm < minElt) & (maxCol > i) ){
        tempIdx <- pivot[maxCol]

        X <- colswap(X, p, i, maxCol)

        pivot[maxCol] = pivot[i]
	pivot[i] = tempIdx

        maxCol <-  maxCol - 1
        colNorm <- ColNorms[i]
        
      }

      if(colNorm >= minElt){
        blockEnd <- blockEnd + 1
      }
      
    }

    rk <-  rk + blockEnd

    for(a in 1:blockSize)
      for(b in 1:rowskp)
        V[b,a] <- 0

    pVcol <- c(1,1)
    pVdiag <- c(1,1)
    pXmdiag <- pXmBlock
    for(colIdx in 1:blockEnd){

      for(a in 1:(rowskp - colIdx))
        V[pVdiag[1] + a - 1, pVdiag[2]] <- X[pXmdiag[1] + a - 1, pXmdiag[2]]

      v1 <- V[pVdiag[1], pVdiag[2]]
      v1abs <- abs(v1)

      if(kp == maxRow)
        qrAux[kp] <- v1abs
      else{
        normV <- 0
        for(a in 1:(rowskp - colIdx)){
          diag <- V[pVdiag[1] + a - 1, pVdiag[2]]
          normV <-  normV + diag*diag
        }

        normV <- sqrt(normV)
        recipNormV <- 1.0 / normV
        qrAux[kp] <- 1.0 + v1abs*recipNormV
        if(v1>=0)
          sgn <- 1.0
        else
          sgn <- -1.0
        scales[colIdx] <- recipNormV*sgn

        fac <- 1.0 + normV/v1abs

        V[pVdiag[1], pVdiag[2]] <- V[pVdiag[1], pVdiag[2]] * fac

        Beta[colIdx] <- -2.9 / (normV * normV + v1abs * (-1.0 + fac * fac));

        Xnew <- X[pXmBlock[1]*(1:rowskp), pXmBlock[2]*(1:min(blockSize, colskp))]
        Vnew <- t(t(V[pVcol[1]*(1:rowskp), pVcol[2]]))
        u <- cXtv(Beta[colIdx], rowskp, min(blockSize, colskp), Xnew, Vnew)

        Xnew <- X[pXmBlock[1]*(1:rowskp), pXmBlock[2]*(1:min(blockSize, colskp))]
        X[pXmBlock[1]*(1:rowskp), pXmBlock[2]*(1:min(blockSize, colskp))] <-
          Xpvut(rowskp, min(blockSize, colskp), Xnew, Vnew, u)
                          
      }

      pVcol <- pVcol + c(0,1)
      pVdiag <- pVdiag + c(1,1)
      pXmdiag <- pXmdiag + c(1,1)
      kp <-  kp + 1
    }

    if (bc < blockCount  && blockEnd > 0) {
      pTcol <- c(1,1)
      pWcol <- c(1,1)
      pVcol <- c(1,1)
      T <- uVtV(rowskp, blockSize, V)

      for(m in 1:blockSize){

        for(a in 1:rowskp){
          W[pWcol[1] + a - 1, pWcol[2]] <- Beta[m]*V[pVcol[1] + a - 1, pVcol[2]]
        }

        if(m>0){
          for(a in 1:rowskp){
            sum <-  0
            for(b in 1:m){
              sum <-  sum + W[b,a] * T[pTcol[1], pTcol[2]+b-1]
            }
            W[pWcol[1], pWcol[2]+a-1] <- W[pWcol[1], pWcol[2]+a-1] + sum
          }
        }
        pWcol <- pWcol + c(0,1)
        pVcol <- pVcol + c(0,1)
        pTcol <- pTcol + c(0,1)
      }

      for(a in 1:blockSize){
        for(b in 1:(colskp - blockSize)){
          sum <- 0
          for(c in 1:rowskp)
            sum <-  sum + W[c,a] *
              X[pXmBlock[1] + c - 1, pXmBlock[2] + blockSize + rowskp - 1]
          WtR[a, b] <- sum
        }
      }

      for(a in 1:rowskp){
        for(b in 1:(colskp - blockSize)){
          sum <- 0
          for(c in 1:blockSize){
            sum <- sum + V[a, c]*WtR[c, b]
          }
          X[pXmBlock[1] + a - 1, pXmBlock[2] + blockSize + rowskp - 1] <-
            X[pXmBlock[1] + a - 1, pXmBlock[2] + blockSize + rowskp - 1] + sum
        }
      }
    }

    pVdiag <- c(1,1)
    pXmdiag <- pXmBlock

    for(l in 1:blockEnd){



      for(a in 1:(rowskp - l - 1)){
        V[pVdiag[1] + 1 + a - 1, pVdiag[2]] <-
          V[pVdiag[1] + 1 + a - 1, pVdiag[2]] * scales[1]
        X[pXmdiag[1] + 1 + a - 1, pXmdiag[2]] <-
          X[pXmdiag[1] + 1 + a - 1, pXmdiag[2]] * V[pVdiag[1] + 1 + a - 1, pVdiag[2]]
      }
      
      pVdiag <-  pVdiag + c(1,1)
      pXmdiag <- pXmdiag + c(1,1)
      
    }

    pXmBlock <- pXmBlock + blockSize * c(1,1)
    colskp <- colskp - blockSize
    rowskp <- rowskp - blockSize
    
  }

  rank <- rk

  return(X)
  
}

n <- 10
p <- 2
X <- matrix(0, ncol=p, nrow=n)
blockSize <- 2^7

for(i in 1:n){
  for(j in 1:p){
    if(j==1){
      X[i,j] <- 1
    }
    else{
      X[i,j] <- i*j
    }
  }
}

Y <-X%*%t(t(c(1:p)))
ColNorms <- getColNorms(n, p, X)

u <- rep(0,n)
T <- matrix(0, blockSize, blockSize)
W <- matrix(0, ncol=blockSize, nrow=n)
V <- W
rank <- 0
WtR <- matrix(0, nrow=blockSize, ncol=-(p-blockSize))
pivot <- rep(0,p)
qrAux <- pivot
Beta <- rep(0, blockSize)
scales <- Beta

R <- GetQRDecompBlocked(n, p, X, ColNorms, blockSize, scales, Beta,
			      qrAux, pivot, WtR, rank, V, W, T, u);




for(i in 1:n)
  for(j in 1:p)
      if(i>j)
        R[i,j] <- 0


R2 <- gpuQr(X)$qr
attr(R2, "Csingle") <- NULL

for(i in 1:n)
  for(j in 1:p)
      if(i>j)
        R2[i,j] <- 0

bR <- solve(t(R)%*%R)%*%t(X)%*%t(t(Y))
bR2 <- solve(t(R2)%*%R2)%*%t(R2)%*%t(t(Y))

lmb <- lm(Y~X[,2])
blm <- lmb$coef
