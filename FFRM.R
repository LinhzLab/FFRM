########################
#### load libraries ####
########################




######    Identifiability condition on U and V (alpha)   #####
ide_U <- function(x) {
  U <- x
  for (i in 1 : nrow(U)) U[i,] <- U[i,] / sqrt( sum(U[i,]^2) )
  svd_U <- svd(U)
  U <- svd_U$u % * % diag(svd_U$d)
  return(U)
}

signal_U <- function(alpha,x) {
  U <- x
  d <- ncol(x)
  for (i in 1:d) if ( sign(alpha[i,1])<0 ) U[i,] <- -U[i,] 
  return(U)
}

ide_V <- function(x) {
  V_qr <- qr( t(x) )
  V <- t( qr.Q( V_qr ) )
  for (i in 1:r) V[i,] <- sign( V[i, (which(V[i,] != 0))[1] ] ) * V[i,]
  return(V)
}



######    function to generate B-spline functions and the derivatives   #####
t_gene <- function(x, break_seq) {
  B <- array(NA,c(n, K1, d))
  Bd <- array(NA,c(n, K1,d))
  for (i in 1 : d) {
    B[,,i] <- bsplineS(x[,i],breaks = break_seq[i,],norder = K1 + 2 - length(break_seq[i,]),nderiv = 0) 
    Bd[,,i] <- bsplineS(x[,i],breaks = break_seq[i,],norder = K1 + 2 - length(break_seq[i,]),nderiv = 1)  
  }
  result <- list(B = B, Bd = Bd)
  return(result)
}



######    A instrument function to help to calculate the derivative   #####
tBF <- function(ta, Ba_d, F_score) {
  tBF <- array(0,c(n, m+K ,d ))
  for (j in 1 : d) {
    for (i in 1 : n) {
      tBF[i,,j] <- c( t( ta[j,] ) %*%  Ba_d[i,,j] ) * F_score[i,]
    }
  }
  return(tBF)
}








FFRM <- function (F_score, Y, n, K, m, d, r, num.spline, u_ini, v_ini, Ba_ini, Ba_d_ini, ta_ini, break_seq, kmax, learning_rate_Ba,learning_rate_u,learning_rate_v, end_condition, lambda = 1e-3, bd = 0.2) {
  

  # F_score : latent score from Factor model and FaFPCA (FaFPCA: https://github.com/LinhzLab/FaFPCA), see feature_extract.R
  # sample.len: sample sizes 
  # K: dimension of latent score from functional predictors
  # m: dimension of latent score from scalar predictors
  # d: number of component funtions
  # r: rank of uv
  # num_spline : number of spline basics
  # u_ini, v_ini : initial value of U,V,alpha
  # Ba_ini, Ba_d_ini : inial values of spline basics functions M_2(alpha[k,] %*%F_score) and the 1-order derivatives
  # ta_ini : initial value of spline coefficeints
  # break_seq : spline knots
  # kmax : maximum number of iteration steps
  # lambda : tuning parameter in adaptive lasso
  # bd : bandwidth of kernel function
  
  K1 <- ncol(ta_ini) # number of spline functions 
  
  psi_ini <- matrix(0, n, d) # set initial value of the estimated mean 
  g_ni <- rep(0, n)
  for (i in 1 : d) {
    psi_ini[,i] <- Ba_ini[,,i] %*% ta_ini[i,]
    g_ini <- g_ini + psi_ini[,i]
  }
  
  U_ini <- matrix(0, n, n) # set the initial of kernel matrix and its 1-order derivatives
  for (i in 1 : n) {
    for (j in 1 : n) {
      U_ini[i,j] <-  (Y[i] - g_ini[i] - Y[j] + g_ini[j]) / bd
    }
  }
  Kern <- 1 / sqrt( 2 * pi ) * exp( - U ^ 2 / 2)
  Kern_d <- 1 / sqrt( 2 * pi ) * exp( - U ^ 2 / 2) * ( -U ) 
  
  
  ###  set initival value  ###
  v <- v_ini
  v <- ide_V(v)
  u <- u_ini
  u <- ide_U(u)
  u <- signal_U(u%*%v,u)
  alpha <- u_ini %*% v_ini
  ta <- ta_ini
  Ba <- Ba_ini
  Ba_d <- Ba_d_ini
  f <- F_score %*% t( alpha_ini )
  psi <- psi_ini
  g <- g_ini
  U <- U_ini
  tBF1 <- tBF(ta, Ba_d, F_score)
  
  
  
  aalpha <- alpha + 1
  uu <- u + 1
  vv <- v + 1
  error_alpha <- array(0, c(d, m+K, kmax+1) )
  error_u <- array(0, c(d, r, kmax+1) )
  error_v <- array(0,c(r, m+K, kmax+1))
  
  
  ###   First step : LASSO estimator   ###
  for (k in 1 : kmax) {  
    ###   Updating ta and calculate the gradient of alpha   ###
    alpha_dev <- matrix(0, d, m+K)
    for (I in 1 : d) {
      a <- rep(0, m+K)
      aa <- rep(0, K1)
      for (i in 1:n) {
        A1 <- sum( Kern[,i] )
        A2 <- rep(0, m+K*q)
        A3 <- rep(0, K1)
        for (j in 1 : n) {
          A2 <- A2 + Kern_d[j,i] * (tBF1[i,,I] - tBF1[j,,I])
          A3 <- A3 + Kern_d[j,i] * (Ba_d[i,,I] - Ba_d[j,,I])
        }
        A2 <- A2 / A1
        A3 <- A3 / A1
        a <- a + A2
        aa <- aa + A3 
      } 
      alpha_dev[I,] <- a / (n*bd)
      ta[I,] <- ta[I,] + 1e-3 *abs( ta[I] ) * aa / ( n*bd )
    }
    
    
    ###   Updating v   ###   
    v_dev <- t(u) %*% alpha_dev
    u_dev <- alpha_dev %*% t(v)
    tt <- 1e-3
    for (l in 1 : (m+K) ) {
      bl1 <- v_dev[,l]
      S <- v[,l] + tt * ( bl1)
      if ( sqrt( sum(S^2) ) <= tt * lambda ) { v[,l] <- rep(0,r) } 
      else { v[,l] <- max( ( 1 - (tt*lambda) / sqrt( sum(S^2) ) ), 0 ) * S}
    }
    v <- ide_V(v)
    
    
    ###   Updating u  ###   
    for (i in 1 : d) {
      u[i,] <- u[i,] + tt * u_dev[i,]
      u[i,] <- u[i,] / sqrt( sum(u[i,]^2) )
    }
    u <- ide_U(u)
    u <- signal_U(u %*% v, u)
    
    
    ###   Other updating   ###
    alpha <- u %*% v
    f <- F_score %*% t(alpha)
    Ba <- t_gene(f,break_seq)$B
    Ba_d <- t_gene(f,break_seq)$Bd
    g <- 0
    for (i in 1 : d) {
      psi[,i] <- Ba[,,i] %*% ta[i,]
      g <- g + psi[,i]
    }
    for (i in 1 : n){
      for (j in 1 : n){
        U[i,j] <-  (Y[i] - g[i] - Y[j] + g[j]) / bd
      }
    } 
    Kern <- 1 / sqrt( 2 * pi ) * exp( - U ^ 2 / 2)
    Kern_d <- 1 / sqrt( 2 * pi ) * exp( - U ^ 2 / 2) * ( -U ) 
    tBF1 <- tBF(ta, Ba_d, F_score)
    error_alpha[,,k] <- abs(aalpha - alpha)
    error_u[,,k] <- norm((uu - u),type = 2)
    error_v[,,k] <- norm((vv - v),type = 2)
    aalpha <- alpha
    uu <- u
    vv <- v
    if (  max(abs(error_alpha[,,k])) <= end_condition)  { break }
    
  }
  
  
  
  ###   Second step : find adaptive LASSO weight  ###
  weight_lasso <- rep(1e8, m+K)
  for (k in 1 : (m+K) ) {
    if ( sqrt( sum(v[,k]^2) )  != 0 ) { weight_lasso[k] <- 1/sqrt( sum(v[,k]^2) ) }
  }
  
  
  ###   Third step : adaptive LASSO estimator   ###
  ###   set initival value  ###
  v <- v_ini
  v <- ide_V(v)
  u <- u_ini
  u <- ide_U(u)
  u <- signal_U(u%*%v,u)
  alpha <- u_ini %*% v_ini
  ta <- ta_ini
  Ba <- Ba_ini
  Ba_d <- Ba_d_ini
  f <- F_score %*% t( alpha_ini )
  psi <- psi_ini
  g <- g_ini
  U <- U_ini
  Kern <- 1 / sqrt( 2 * pi ) * exp( - U ^ 2 / 2)
  Kern_d <- 1 / sqrt( 2 * pi ) * exp( - U ^ 2 / 2) * ( -U ) 
  tBF1 <- tBF(ta, Ba_d, F_score)
  aalpha <- alpha + 1
  uu <- u + 1
  vv <- v + 1
  error_alpha <- array(0, c(d, m+K, kmax+1) )
  error_u <- array(0, c(d, r, kmax+1) )
  error_v <- array(0,c(r, m+K, kmax+1))
  k_stop <- 0
  
  for (k in 1 : kmax) {
    ###   Updating ta and calculate the gradient of alpha   ###
    alpha_dev <- matrix(0, d, m+K)
    for (I in 1 : d) {
      a <- rep(0, m+K)
      aa <- rep(0, K1)
      for (i in 1 : n) {
        A1 <- sum( Kern[,i] )
        A2 <- rep(0, m+K)
        A3 <- rep(0, K1)
        for (j in 1 : n) {
          A2 <- A2 + Kern_d[j,i] * (tBF1[i,,I] - tBF1[j,,I])
          A3 <- A3 + Kern_d[j,i] * (Ba_d[i,,I] - Ba_d[j,,I])
        }
        A2 <- A2 / A1
        A3 <- A3 / A1
        a <- a + A2
        aa <- aa + A3 
      } 
      alpha_dev[I,] <- a / (n*bd)
      ta[I,] <- ta[I,] + learning_rate_Ba[I] * abs(ta[I,]) * aa/(n*bd)
    }
    
    
    
    ###   Updating v  ###  
    u_dev <- alpha_dev %*% t(v)
    v_dev <- t(u) %*% alpha_dev
    for (l in 1 : (m+K) ) {
      bl1 <- v_dev[,l]
      S <- v[,l] + learning_rate_v[l] * ( bl1)
      if ( sqrt(sum(S^2)) <= learning_rate_v[l] * lambda * weight_lasso[l] ) { v[,l] <- rep(0,r) } 
      else { v[,l] <- max( ( 1-(learning_rate_v[l] * lambda * weight_lasso[l]) / sqrt( sum(S^2) ) ), 0 ) * S}
    }
    v <- ide_V(v)
    
    
    
    ###   Updating u  ###  
    for (i in 1:d) {
      u[i,] <- u[i,] + learning_rate_u * u_dev[i,]
      u[i,] <- u[i,]/sqrt( sum(u[i,]^2) )
    }
    u <- ide_U(u)
    u <- signal_U(u%*%v, u)
    
    
    
    ###   Other updating   ###
    alpha <- u %*% v
    f <- F_score %*% t(alpha)
    Ba <- t_gene(f, break_seq)$B
    Ba_d <- t_gene(f, break_seq)$Bd
    g <- 0
    for (i in 1 : d) {
      psi[,i] <- Ba[,,i] %*% ta[i,]
      g <- g + psi[,i]
    }
    for (i in 1 : n){
      for (j in 1 : n){
        U[i,j] <-  (Y[i] - g[i] - Y[j] + g[j]) / bd
      }
    } 
    Kern <- 1 / sqrt( 2 * pi ) * exp( - U ^ 2 / 2)
    Kern_d <- 1 / sqrt( 2 * pi ) * exp( - U ^ 2 / 2) * ( -U ) 
    tBF1 <- tBF(ta, Ba_d, F_score)
    error_alpha[,,k] <- abs(aalpha - alpha)
    error_u[,,k] <- norm((uu - u),type = 2)
    error_v[,,k] <- norm((vv - v),type = 2)
    aalpha <- alpha
    uu <- u
    vv <- v
    if (  max( max(abs(error_u[,,k])), max(abs(error_v[,,k])) ) <= end_condition )  {break}  
  }
  
  
  
  ###   output result   ###
  return( list(u, v, alpha, ta, Ba, Ba_d, psi, g) ) 
  
  
  
}