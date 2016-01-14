gmm_s <- function(y, X, Z){
  # _s indicates the structural equation
  # use "lars" to compute the path of the lasso solutions
  
  D <- dim(X)
  N <- D[1]
  K <- D[2] 

  PZ <- solve( t(Z) %*% Z )
  svd.PZ <- svd(PZ)
  HProj <- diag(sqrt(svd.PZ$d)) %*% t(svd.PZ$v) %*% t(Z)

  PX = HProj %*% X
  Py = HProj %*% y

  reg <- lars(PX, Py, type = "lasso", normalize = F, intercept = F)


  bhat <- coef(reg, newx = X, s = seq(0.01,1,length.out = 100), mode = "fraction")
  bhat <- t(bhat)

  return(bhat) 
}

##############################

two.equation.est <- function(k0, y, X, Z){
  
  bhat.s <- gmm_s(y, X, Z)
  bhat.r <- gmm_s(X[,k0], X[, -k0], Z)
  

  #### automated ######
  mA.s <- mIC(y, X, Z, bhat.s, "AIC")
  choice.s <- which.min(mA.s)
  bhat.ss <- bhat.s[,  choice.s]
  
  mA.r <- mIC(X[,k0], X[,-k0], Z, bhat.r, "AIC")
  choice.r <- which.min(mA.r)
  bhat.rr <- bhat.r[, choice.r ]
  ###############
  
  bhat <- list(s = bhat.ss, r = bhat.rr)
}
########################################

Proj <-  function(Z) { 
  y <- Z %*% solve( t(Z) %*% Z ) %*% t(Z)
  return(y)
}
##############################
stat <- function(k0, y, X, Z, bb, display = F){
  # bb: is a outcome of "two.equations"
  
  # from the reduced-from

  Omega <- Proj(Z)
  v <- X[,k0] - X[,-k0] %*% bb$r # the residual
  vPX <- t(v) %*% Omega %*% X[,k0]; # v * P * X1
        
  # the part for the standard error of the estimate
  tau <- sqrt(t(v) %*% Omega %*% v) / vPX; 
  
  # from the structural form
  e <- as.vector( y - X %*% bb$s )
  sig.e <- sd(e)
  if(display == T ) { cat("outside sig.e = ", sig.e, "\n") }
  #######################
    b.o = bb$s[k0] # the original lasso estimate
    corr = t(v) %*% Omega %*% e / vPX # the bias correction

    b.corr <- b.o + corr # the corrected beta estimate
          
    CI.L = b.corr - 1.96 * sig.e * tau # the lower end of CI
    CI.U = b.corr + 1.96 * sig.e * tau # the upper end of CI
    CI <- c(CI.L, CI.U) # confidence interval
  return( c(b0 = b.o, b.corr = b.corr, CI = CI, sig.e = sig.e) )
  }

##############
mIC <- function(y, X, Z, bhat, method){
  # the modified information criterion to choose the tuning parameter
  
  N <- dim(X)[1]
  K <- dim(X)[2]
  C <- log(log(K))
  
  y <- as.vector(y)
  
  PZ <- Z %*% solve( t(Z) %*% Z ) %*% t(Z)
  e <- y - X%*%bhat
  sigmaS <- diag(1/N * t(e) %*% PZ %*% e )  
  ###############
  S <- colSums( abs(bhat) > 0.0001 )
  if (method == "BIC") {
    
    crit <- sigmaS + S * log(N)*C/N
  } else if (method == "AIC") {
    crit <- sigmaS + S * 2 * C/N
  }
  return(crit)
}

#######################
LIML <- function(y, X, Z, method){
  # method == either "tsls" or "liml"
  
  # when exactly identied, LIML is exact the same as 2SLS
  # TSLS and LIML are different when overidentied.
  N <- dim(Z)[1]
  
  I <- diag(rep(1,N))
  PZ <- Z %*% solve( t(Z) %*% Z ) %*% t(Z)
  M <- I - PZ
  
  # to determine "k" of the k-class estimator
  if (method == "liml") {
    # I follow the closed form in Hayashi
    bigY <- cbind(y, X)
    W.tilde <- crossprod(bigY)
    W <- t(bigY) %*% M %*% bigY
    
    kk <- eigen(W.tilde %*% solve(W), only.values = T)$values
    k <- min( abs(kk) )
   
 #    print( k )
  } else if (method == "tsls"){
    k <- 1
  }
  
  beta <- solve( t(X) %*% (I - k * M) %*% X ) %*% (t(X) %*% (I - k * M) %*% y )
  return(beta)
}
