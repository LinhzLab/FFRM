###   Exampe 1  :  FFRM   ###

####  The generation of U,V   ####
v.true <- matrix(0, r, m+K*q)
v.true[,c(1, 2, m+1, m+2)] <- c(1,1)
v.true[,c(3, 4, m+3, m+4)] <- c(1,-1)
for (i in 1 : r) v.true[i,] <- v.true[i,] / sqrt( sum(v.true[i,] ^ 2) )
u.true <- matrix(0, d, r)
u.true[,1] <- c(3, 3, 2, 2)
u.true[,2] <- c(2, -2, -1, 1)
for (i in 1 : d) u.true[i,] <- u.true[i,] / sqrt( sum(u.true[i,] ^ 2) )
u.true <- svd(u.true)$u %*% diag(svd(u.true)$d)
u.true <- signal.U(u.true%*%v.true,u.true)

####  The generation of f=alpha.1%*%Fi+alpha.2%*%score.i ####
alpha.true <- u.true %*% v.true
f.true <- F.score.true %*% t(alpha.true) 

###   The generation of psi(f.true)   ####
psi.true1 <- abs( f.true[,1] + 1 )
psi.true2 <- ( f.true[,2] ) ^ 2
psi.true3 <- sin( f.true[,3] )
psi.true4 <- cos( f.true[,4] )
g.true <- psi.true1 + psi.true2 + psi.true3 + psi.true4

###   The generation of Y   ###
Y.error <- rnorm(n = n, 0, 0.2) # setting 1

sel <- rbinom(n, size = 1, prob = 1/3) 
Y.error = sel * rnorm(n, -4, 2/5) + (1-sel) * rnorm(n, 2, 1/5) # setting 2


Y.error <- abs( f.true[,1] +1) * rnorm(n = n, 0, 0.1) # setting 3
Y <-  g.true + Y.error



###   Exampe 2   ###

g.true <- rep(0,n)
for (j in 1 : 10) {
  g.true <- g.true + cos( F.score.true[,j] )
  g.test.true <- g.test.true + cos( F.score.test.true[,j] )
}
for (j in 11 : 20) {
  g.true <- g.true + sin(F.score.true[,j])
  g.test.true <- g.test.true + sin( F.score.test.true[,j] )
}
Y.error <- rnorm(n,0,0.5)
Y  <- g.true + Y.error



###   Exampe 3  :  kong2016   ###
###   generation true uFPCA score   ###
score.new.true <- array(0,c(n, K*q, p) )
for (j in 1 : p) {
  bjq <- rep(B.true[j,], each = K)
  for (i in 1 : n) {
    score.new.true[i,,j] <- bjq * score.true[i,]
  }
}

###   generation Y   ###
gamma.true <- rep(0, M)
gamma.true[1 : 20] <- 0.1/sqrt(10)
beta.true <- matrix(0, p, K*q)
for (j in 1:20) {
  beta.true[j,] <- 0.1
}
g.true <- rep(0,n)
for (i in 1:n) {
  for (j in 1:p) {
    g.true[i] <- g.true[i] + sum(beta.true[j,] * score.new.true[i,,j])
  }
  g.true[i] <- g.true[i] + sum(gamma.true * Z[i,])
}
Y  <- g.true + rnorm(n, 0, 0.5)


###   Exampe 4  :  wong2019   ###
###   generation true mFPCA score and change ie to zeta to generate Y   ###
ufpca.score.combine.true <- matrix(0, n, p*K*q)
for (i in 1 : n ) {
  for (j in 1:  p) {
    ufpca.score.combine.true[i, ( (j-1) * (K*q) + 1) : ( (j-1) * (K*q) +K * q) ] <- score.new.true[i,,j]
  }
}
m.svd.true <- svd(t(ufpca.score.combine.true) %*% (ufpca.score.combine.true))
mfpca.score.true <- matrix(0, n, 10)
for (i in 1 : n) mfpca.score.true[i,] <- m.svd.true$u[1 : 10,] %*% ufpca.score.combine.true[i,] 

var.mfpca.score.true <- var(rbind(mfpca.score.true,mfpca.score.test.true))
lambda.mfpca.score.true <- svd(var.mfpca.score.true)$d
mfpca.zeta.true <- matrix(0, n, 10)
for (k in 1 : 10) {
  for (i in 1 : n) {
    mfpca.zeta.true[i,k] <- pnorm(1 / sqrt( var(mfpca.score.true[,k]) ) * mfpca.score.true[i,k], 0, 1)
  }
}

###   generation Y   ###
f1.true <- 3 * mfpca.zeta.true[,1] - 3/2
f2.true <- sin( 2 * pi * (mfpca.zeta.true[,3] - 0.5) )
f3.true <- 4 * ( mfpca.zeta.true[,5] - 1/2 ) ^ 3 - 8/9
f4.true <- 2 * ( cos(mfpca.zeta.true[,7] + 0.5) )
theta.Z.true <- rep(0, p)
theta.Z.true[1 : 10] <- c(0.1, 0.1, 0.1, 0.1, 0.1, -0.1, -0.1, -0.1, -0.1, -0.1) / sqrt(10)
g.true <- Z %*% theta.Z.true
g.true <- f1.true + f2.true + f3.true + f4.true + g.true
Y  <- g.true + rnorm(n, 0, 0.5)