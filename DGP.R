# only contains functions that generates the data

DGP1 <- function(N, K, rho, K0, K.red, b1, b2, p1, p2){
  
  b0 <- rep(b2, K)
  b0[1:K0] <- b1
  
  E <- matrix( rnorm( (K+1)*N ), ncol = N)
  U <- as.vector( E[1,] )
  V <- t( E[ 2:(K+1), ] )
  
  
  Z1 <- matrix( rnorm(N*K, 0, 1), nrow = N )
  Z2 <- matrix( rnorm(round( N*K * (K.red - 1) ), 0, 1), nrow = N )
  # Z3 <- matrix( rnorm(N*K, 0, 0.6), nrow = N )
  Z1 <- normalize( Z1 )
  z2 <- normalize( Z2 )
  # browser()
  X <- p1 * Z1  + rho * U + p2 * V
  X <- normalize(X)
  
  y = X %*% b0 + U;
  return( list( y = y , X = X, Z = cbind(Z1, Z2), b0 = b0, e0 = U ) )
}

###############################

DGP2 <- function(N, K, rho, K0, K.red, b1, b2, p1, p2){
  
  b0 <- rep(b2, K)
  b0[1:K0] <- b1
  
  E <- matrix( rnorm( (K+1)*N ), ncol = N)
  U <- as.vector( E[1,] )
  V <- t( E[ 2:(K+1), ] )
  
  
  Z1 <- matrix( rnorm(N*K, 0, 1), nrow = N )
  Z2 <- matrix( rnorm(round( N*K * (K.red - 1) ), 0, 1), nrow = N )
  # Z3 <- matrix( rnorm(N*K, 0, 0.6), nrow = N )
  Z1 <- normalize( Z1 )
  z2 <- normalize( Z2 )
  
  U1 <- matrix( rep(rho * U, K0), ncol = K0 )
  U2 <- matrix( rho * rnorm(N * (K-K0)), ncol = K - K0)
  
  X <- p1 * Z1  + p2 * V +  cbind(U1, U2)
  X <- normalize(X)
  
  y = X %*% b0 + U;
  return( list( y = y , X = X, Z = cbind(Z1, Z2), b0 = b0, e0 = U ) )
}

###############################


DGP3 <- function(N, K, rho, K0, K.red, b1, b2, p1, p2){
  
  b0 <- rep(b2, K)
  b0[1:K0] <- b1
  
  E <- matrix( rnorm( (K+1)*N ), ncol = N)
  U <- as.vector( E[1,] )
  V <- t( E[ 2:(K+1), ] )
  
  
  Z1 <- matrix( rnorm(N*K, 0, 1), nrow = N )
  Z2 <- matrix( rnorm(round( N*K * (K.red - 1) ), 0, 1), nrow = N )
  # Z3 <- matrix( rnorm(N*K, 0, 0.6), nrow = N )
  Z1 <- normalize( Z1 )
  z2 <- normalize( Z2 )
  

  AA <- matrix( 0, ncol = K, nrow = K)
  mult <- sqrt(1/3)
  diag(AA) <- mult * p1
  AA[1, K] <- mult * p1
  AA[K, 1] <- mult * p1
  AA[seq(from = 2, to = K*K, by = K+1)] <- mult * p1
  AA[seq(from = K+1, to = K*K, by = K+1)] <- mult * p1
  # browser()
  U1 <- matrix( rep(rho * U, K0), ncol = K0 )
  U2 <- matrix( rho * rnorm(N * (K-K0)), ncol = K - K0)
  
  X <-  Z1%*%AA  + p2 * V +  cbind(U1, U2)
  X <- normalize(X)
  
  y = X %*% b0 + U;
  return( list( y = y , X = X, Z = cbind(Z1, Z2), b0 = b0, e0 = U, AA = AA ) )
}

###############################

DGP4 <- function(N, K, rho, K0, K.red, b1, b2, p1, p2, AA = diag(K)){
  
  b0 <- rep(b2, K)
  b0[1:K0] <- b1
  
  E <- matrix( rnorm( (K+1)*N ), ncol = N)
  U <- as.vector( E[1,] )
  V <- t( E[ 2:(K+1), ] )
  
  
  Z1 <- matrix( rnorm(N*K, 0, 1), nrow = N )
  Z2 <- matrix( rnorm(round( N*K * (K.red - 1) ), 0, 1), nrow = N )
  # Z3 <- matrix( rnorm(N*K, 0, 0.6), nrow = N )
  Z1 <- normalize( Z1 )
  z2 <- normalize( Z2 )
  
  A <-  sqrt( 0.5^AA ) 

  U1 <- matrix( rep(rho * U, K0), ncol = K0 )
  U2 <- matrix( rho * rnorm(N * (K-K0)), ncol = K - K0)
  
  X <- Z1 * p1 +  V%*%A +  rho * U #  cbind(U1, U2)
  X <- normalize(X)
  
  y = X %*% b0 + U;
  return( list( y = y , X = X, Z = cbind(Z1, Z2), b0 = b0, e0 = U, AA = AA ) )
}

###############################

DGP6 <- function(N, K, rho, K0){
  
  b0 <- c( rep(1, K0),  rep(0.00, K-K0) )
  
  
  E <- matrix( rnorm( (K+1)*N ), ncol = N)
  U <- as.vector( E[1,] )
  V <- t( E[ 2:(K+1), ] )
  
  Z <- matrix( rnorm(N*K, 0, 0.6), nrow = N )
  Z <- normalize( Z )
  X <- 1 + 1 * Z + rho * U + 0.1 * V
  X <- normalize(X)
  
  y = X %*% b0 + U;
  return( list( y = y , X = X, Z = Z, b0 = b0, e0 = U ) )
  }



DGP8 <- function(N, K, rho, K0, b1, b2){
# works well with the following parameters
#   #############################
#   N <- 160
#   K <- 40
#   K0 <- 10
#   rho <- 0.3
#   
#   b1 <- 1
#   b2 <- 0.0
#   
#   R <- 20
# 5-fold CV
  
  b0 <- c( rep(b1, K0),  rep(b2, K-K0) )
  
  
  E <- matrix( rnorm( (K+1)*N ), ncol = N)
  U <- as.vector( E[1,] )
  V <- t( E[ 2:(K+1), ] )
  
  K.red <- K * 1
  
  Z1 <- matrix( rnorm(N*K, 0, 0.6), nrow = N )
  Z2 <- matrix( rnorm(N*K.red, 0, 0.6), nrow = N )
  # Z3 <- matrix( rnorm(N*K, 0, 0.6), nrow = N )
  Z1 <- normalize( Z1 )
  z2 <- normalize( Z2 )
  # browser()
  X <- 1 + 0.5 * (Z1 + Z2) + rho * U + 0.1 * V
  X <- normalize(X)
  
  y = X %*% b0 + U;
  return( list( y = y , X = X, Z = cbind(Z1, Z2), b0 = b0, e0 = U ) )
}


normalize<- function(W) {
  nrow <- dim(W)[1]
  W.mean <- colMeans(W)
  W.demean <- W - W.mean
  W.std    <- apply(W, 2, sd)
  X <- W.demean/W.std
  return(X)
}

#############################

normalize<- function(W) {
  nrow <- dim(W)[1]
  W.mean <- colMeans(W)
  W.demean <- W - W.mean
  W.std    <- apply(W, 2, sd)
  X <- W.demean/W.std
  return(X)
}
