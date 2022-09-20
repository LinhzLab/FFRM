########################
#### load libraries ####
########################
library("MASS")
library("fda")





sim.XZ <- function (sample.len, time.len, p, q, K, M, m, F.score.var) {
  
  # sample.len: sample size
  # time.len: number of time points for each individual
  # p: dimension of functional predictors 
  # q: dimension of loadings B in functional predictors 
  # K: number of eigenfunctions 
  # M: dimenson of scalar predictors
  # m: dimension of loadings L in scalar predictors 
  # F.score.var: covariance matrix of F and score 
  
  # generate time points
  
  time <- list()
  
  for (i in 1 : sample.len) {
    set.seed(1)
    time[[i]] <- sort( runif(time.len[i], min = 0, max = 10) )
  }
  
  
  ### generate loading B   ###
  
  ideqr.B <- function(X) {
    B.qr <- qr(X)
    B <- sqrt(p) * qr.Q(B.qr)
    for (i in 1 : q) B[,i] <- sign( B[1,i] ) * B[,i]
    return(B)
  }
  set.seed(5)
  b.var <- matrix(0, p, p)
  for (i in 1 : p) {
    for (j in 1 : p) {
      b.var[i,j] <- 0.5 ^ ( abs(i - j)) 
    }
  }
  b.ran <- mvrnorm(sample.len, rep(0,p), b.var)
  Bn.true <- b.ran %*% t(b.ran) 
  Bn.qr <- qr(Bn.true)
  Bn.true <- t(b.ran) %*% qr.Q(Bn.qr)[,1 : q] 
  B.true <- sqrt(p) * qr.Q( qr(Bn.true) )
  B.true  <- ideqr.B(B.true)
  
  
  ###   generate eigenfunctions ###
  
  eigenfun <- list()
  for (i in 1 : sample.len) {
    eigenfun[[i]] <- array(0,c(time.len[i], K, q))
    
    for (l in 1 : q)  {
      for (k in 1 : (K/2)) {
        eigenfun[[i]][,(2 * k) - 1,l] <- sin( (2 * k - 1) * pi * time[[i]] / 10) / sqrt(5)
        eigenfun[[i]][,2 * k,l] <- cos((2 * k - 1)  * pi * time[[i]] / 10) / sqrt(5)
      }
    }
    

  }
  
  
  ### generate loading L   ###
  
  idesvd.L <- function(x) {
    L.svd <- svd(x)
    L <- L.svd$u %*% diag( L.svd$d )
    for (i in 1 : m) L[,i] <- sign( L[1,i] + 0.1 ) * L[,i]
    return(L)
  }
  set.seed(5)
  l.var <- matrix(0, M, M)
  for (i in 1 : M) {
    for (j in 1 : M) {
      l.var[i,j] <- 0.5 ^ ( abs(i - j) )
    }
  }
  l.ran <- mvrnorm(sample.len, rep(0,M), l.var)
  Ln.true <- l.ran %*% t(l.ran) 
  Ln.qr <- qr(Ln.true)
  Ln.true <- t(l.ran) %*% qr.Q(Ln.qr)[, 1 : m] 
  Ln.svd <- svd(Ln.true)
  L.true <- Ln.svd$u %*% diag(Ln.svd$d)
  b.ran <- mvrnorm(sample.len, rep(0,p), b.var)
  Bn.true <- b.ran %*% t(b.ran) 
  Bn.qr <- qr(Bn.true)
  Bn.true <- t(b.ran) %*% qr.Q(Bn.qr)[,1 : q] 
  B.true <- sqrt(p) * qr.Q( qr(Bn.true) )
  B.true  <- ideqr.B(B.true)
  
  
  
  ####   generate score and F   ###
  
  idesvd.score <- function(x) {
    score <- matrix(0, sample.len, q * K)
    score.n <- array(NA, c(q, K, sample.len))
    score.n <- x
    for (i in 1 : sample.len) score[i,] <- as.vector( t(score.n[,,i]) )
    score.svd <- svd(score)
    score <- score.svd$u %*% diag(score.svd$d)
    for (j in 1 : ( q * K )) score[,j] <- sign( score[1,j] ) * score[,j]
    for (i in 1 : sample.len) score.n[,,i] <- t(matrix(score[i,], K, q))
    result <- list(score.n = score.n, score = score)
    return(result) 
  }
  

  set.seed(3)
  F.score.true <- mvrnorm(n = sample.len,rep(0, m+K*q), F.score.var)
  F.true.ran <- F.score.true[,1 : m]
  F.true <- ideqr.F(F.true.ran)
  score.true.ran <- F.score.true[,(m+1):(m+K*q)]
  score.n.true.ran <- array(NA,c(q,K,sample.len))
  for (i in 1 : sample.len) score.n.true.ran[,,i] <- t(matrix(score.true.ran[i,],K,q))
  score.n.true <- idesvd.score(score.n.true.ran)$score.n
  score.true <- idesvd.score(score.n.true.ran)$score
  F.score.true <- cbind(F.true,score.true) 
  
  
  
  
  ### generate functional data  ###
  
  Ht.true <- list()
  for (i in 1 : sample.len) {
    
    Ht.true[[i]] <- matrix(0, q, time.len[i])
    for (j in 1 : q) {
      for (t in 1 : time.len[i]) {
        Ht.true[[i]][j,t] <- t( score.n.true[j,,i] ) %*% eigenfun[[i]][t,,j]
      }
    }
    
  }
  
  X <- list()
  
  for (i in 1 : sample.len) X[[i]] <- B.true %*% Ht.true[[i]] + matrix( mvrnorm(n = time.len[i], rep(0, p), diag(0.01, p, p)), 
                                                                        p, time.len[i], byrow = F ) 
  
  ### generate scalar data   ###
  Z <- F.true %*% t(L.true) + matrix(mvrnorm(n = sample.len,rep(0,M), diag(0.01,M,M)), sample.len, M, byrow = T)
  
  
  
  return( list( Z = Z, X = X, B = B.true, score = score.true, eigenfun = eigenfun, L = L.true, F.score = F.score.true) )
  
}




###################################
#### generation of the example ####
###################################

sample.len <- 100 # 500
time.len <- rep(50, 100) 
p <- M <- 100 # 500 
q <- 5
K <- 2
m <- 5

F.score.var <- matrix(0,m + K * q, m + K * q)
for (i in (m + 1) : (m + K * q)) F.score.var[i, i] <- 1/ (i-m)
for (i in 1 : m) {
  for (j in 1 : m) {
    F.score.var[i,j] <- 0.5^( abs(i-j) )
  }
}
F.score.corr <- matrix(0, m, K*q)
for (i in 1 : m) {
  for (j in 1 : ( K * q )) {
    F.score.corr[i,j] <- 0.2 ^ ( abs(i-j) + 1 )
  }
}
F.score.var[1 : m, (m+1) : (m+K*q)] <- F.score.corr
F.score.var[(m+1) : (m+K*q), 1 : m] <- t(F.score.corr)


run.sim <- sim.XZ(sample.len, time.len, p, q, K, M, m, F.score.var)
X <- run.sim$X
Z <- run.sim$Z
F.score <- run.sim$F.score



