# the functions in this script
# gmm_s: estimation using LARS
# gmm_s_pen: estimation using given lambda. "penalized".
# get.lambda: using the "Haste and Park" to obtain the path of lambdas.
# 


gmm_s <- function(y, X, Z){
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
###########################################
gmm_s_pen <- function(y, X, Z, lambda){

  
  D <- dim(X)
  N <- D[1]
  K <- D[2] 
  
  PZ <- solve( t(Z) %*% Z )
  svd.PZ <- svd(PZ)
  HProj <- diag(sqrt(svd.PZ$d)) %*% t(svd.PZ$v) %*% t(Z)
  
  PX = HProj %*% X
  Py = HProj %*% y
  
  reg1 <- penalized( response = Py, penalized =  PX ,  unpenalized = ~0, 
                     maxiter = 20,
                     lambda1 = lambda, model = "linear", trace = F
                     )
  return(reg1)
}
#######################################
get.lambda <- function(y, X, Z, lengt){
  
  D <- dim(X)
  N <- D[1]
  K <- D[2] 

  PZ <- solve( t(Z) %*% Z )
  svd.PZ <- svd(PZ)
  HProj <- diag(sqrt(svd.PZ$d)) %*% t(svd.PZ$v) %*% t(Z)

  PX = HProj %*% X
  Py = HProj %*% y

  # reg <- lars(PX, Py, type = "lasso", normalize = F, intercept = F)

  reg <- penalized( response = Py, penalized =  PX , 
                    unpenalized = ~0, lambda1 = 0.01, model = "linear", 
                    steps = "Park", trace = F )
  
  lambda0 <- sapply(reg, lambda1) # the lambda from "Park"
  e <- sapply(reg, residuals) # get the residual
  multiple <- sqrt( colSums(e^2)/N ) 
  # using the residual to adjust lambda
  # for the purpose of iterative estimation
  lambda.round2 <- multiple * lambda0
  
  # extract a subset of the lambdas. "Park" gives too many points.
  # this is only for the sake of speed in simulation
  lambda0.length <- length(lambda0)
  lambda.index <- seq(from  = 1, to = lambda0.length, length.out = lengt)
  lambda.index <- round(lambda.index)
  lambda.round2 <- lambda.round2[lambda.index]
  
  return(lambda.round2)
}
##############################
two.equation.est <- function(k0, y, X, Z){
  #not useful
  
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
Proj <-  function(Z) { y <- Z %*% solve( t(Z) %*% Z ) %*% t(Z) }
##############################
stat <- function(k0, y, X, Z, bb, display = F){
  # calculate the bias and variance
  # to construct the confidence interval
  
# input:  
  # bb: two vectors. One is the beta of the structural equation. 
  #     the other is the beta of the reduced-form equation.
  
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
  return( c(b0 = b.o, b.corr = b.corr, CI = CI, sig.e = sig.e, tau = tau) )
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
