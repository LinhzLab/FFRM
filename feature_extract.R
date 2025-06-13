########################
#### load libraries ####
########################


Bt.spline <- function (sample.len, time, spline.break, num.spline) {
  Bt <- list()
  for (i in 1 : sample.len) Bt[[i]] <- bsplineS(time[[i]], breaks = spline.break, 
                                                norder = num.spline + 2 - length(spline.break), nderiv = 0)
  for (k in 1 : num.spline) {
    btkl <- c( Bt[[1]][,k] )
    for (i in 2 : sample.len) btkl <- c( btkl, Bt[[i]][,k] )
    for (i in 1 : sample.len) Bt[[i]][,k] <- Bt[[i]][,k] / norm(btkl, type="2") * sqrt( sample.len )
  }
  
  return( Bt = Bt )
}


FaFPCA <- function (X, time, sample.len, time.len, p, q, K, num.spline, Bt, lambda = 1e-3) {
  
  # X : functional predictors
  # time : list from simulation setiings
  # sample.len: sample sizes 
  # time.len: number of time points for each individual
  # p: dimension of predictors
  # q: dimension of loadings
  # K: number of eigenfunctions
  # num.spline : number of spline basics
  # lambda : tuning parameter to guarantee  invertibility 
  
  # estimate loading 
  
  XX <- matrix(0, p, p)
  for (i in 1 : sample.len) XX <- XX + ( X[[i]] %*% t( X[[i]] ) ) / time.len[i]
  XX.svd <- svd(XX)
  B <- XX.svd $ u[,1:q] * sqrt(p)
  for (k in 1 : q) B[,k] <- sign( B[1,k] ) * B[,k]
  
  
  # estimate score and eigenfunction
  
  Theta.q <- array(NA, c(K, num.spline, q)) 
  Theta <- matrix(0, q, K * num.spline)
  score.n <- array(NA, c(q, K, sample.len))
  score <- matrix(0, sample.len, q * K)
  
  for (k in 1 : q) {
    WK <- matrix(0, num.spline, sample.len)
    
    for (i in 1 : sample.len) {
      WK.1 <- solve( t( Bt[[i]] ) %*% Bt[[i]] + lambda * diag( num.spline ) ) 
      WK.2 <- rep(0, num.spline)
      for (t in 1 : time.len[i]) WK.2 <- WK.2 + Bt[[i]][t,] * c( t( B[,k] ) %*% X[[i]][,t] ) / p
      WK[, i] <-  WK.1 %*% WK.2
    }
    
      WK.svd <- svd( WK %*% t(WK) )
      Theta.q[,,k] <- t( WK.svd$u[,1 : K] )
      for (i in 1 : sample.len) score.n[k,,i] <-  c( Theta.q[,,k] %*% WK[,i] )

  }
  #for (k in 1 : (K * q)) score[,k] <- sign( score[1,k] ) * score[,k]
  
  for (i in 1 : sample.len) score[i,] <- as.vector( t( score.n[,,i] ) )
  for (k in 1 : q) Theta[k,] <- as.vector( Theta.q[,,k] )
  
  return( list(B = B, score = score, Theta.q = Theta.q, Theta = Theta) )
  
}



FM <- function (Z, sample.len, M, m) {
  
  # Z : scalar predictors
  # sample.len: sample sizes 
  # M : dimension of scalar predictors
  # m: dimension of factors
  
  Z_eigen <- svd(Z %*% t(Z) / M)
  F <- sqrt(sample.len) * Z_eigen$u[, 1 : m] 
  for (l in 1 : m) F[ , l] <- sign(F[ 1, l]) * F[ , l]
  L <- t(F) %*% Z / sample.len
  return( list(F = F, L = L) )
  
}

score <- FaFPCA(X, time, sample.len, time.len, p, q, K, num.spline, Bt, lambda = 1e-3)
F <- FM(Z, sample.len, M, m)
F_score <- cbind(F, score)