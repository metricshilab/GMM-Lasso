# this script does the hard work of estimation


tsls.transform <- function(y, X, Z){

  PZ <- solve( t(Z) %*% Z )
  svd.PZ <- svd(PZ)
  HProj <- diag(sqrt(svd.PZ$d)) %*% t(svd.PZ$v) %*% t(Z)
  # HProj is a half of the projection matrix Z %*% inv(Z'Z) %*% Z'
  # here using the singular decomposition is correct.
  # instead, it would go wrong if we use the eigenvalue and eigenvector decomposition
  # because the eigen-decomposition will have 
  # A = t(P) %*% D %*% P, but t(P) != P. 
  # eigen decomposition is NOT the half-half solution.  
  
  PX = HProj %*% X
  Py = HProj %*% y

  return( list(PX = PX, Py = Py) )
}



gmm_s <- function(y, X, Z){
    # _s indicates the structural equation
    # use "lars" to compute the path of the lasso solutions
    
    D <- dim(X)
    N <- D[1]
    K <- D[2] 
    

    PXy <- tsls.transform(y, X, Z)
    PX  <- PXy$PX
    Py  <- PXy$Py
    
    PX1 = PX
    Py1 = Py

    reg <- lars(PX, Py, type = "lasso", normalize = F, intercept = F)
    
    bhat <- coef(reg, s = seq(0.01,1,length.out = 100), mode = "fraction")
    bhat <- t(bhat)
    
    
    cvlars <- cv.lars.2sls(x = X, y = y, z = Z, K = 10, plot.it = F, se = F,
                      index = seq(0.01,1,length.out = 100),
                      mode = "fraction", type = "lasso")
    
    b.cv <- coef(reg, s = seq(0.01,1,length.out = 100), mode = "fraction")
    b.cv <- t(b.cv)
    bhat.cv <- b.cv[, which.min(cvlars$cv)]
    
    # browser()
    return( list( bhat = bhat, bhat.cv = bhat.cv)  )
}


############################
summ.stat <- function(B, b0){
    # it is a function called in "master1"
    # calculate the summary statistics for MSE, squared-bias, and variance
    
    biasS <- colSums( (apply( B, c(1,2), mean) - b0)^2 )
    vars <- colSums( apply(B, c(1,2), var) )
    MSE  <- biasS + vars
    tab <-  rbind( biasS, vars, MSE) 
    colnames(tab) <- c("AIC", "BIC", "CV", "GMM", "Infeasible")
    return(tab )
}


######################
dist <- function(K){
    # invoked in "master" 1 about DGP 4
    
    A <- matrix(0, ncol = K, nrow = K)
    
    for ( kkk in 1:K ){
        if (kkk == 1) {A[kkk, ] <- 1:K } 
        else if (kkk == K) {
            A[kkk,  ] <- rev(1:K ) 
        }   else {
            A[ kkk, ] <- c( rev( 1:(kkk-1) )+1 , 1:(K-kkk + 1)  )
        }  
    }
    return(A)
}

#######################################
master1 <- function(N, K, K0, K.red, rho, b1, b2, p1, p2, R, DGP){

    B <- array(0, dim = c(K, 5, R) )
    
    # why here I only have DGP 2, 3 and 4? Where is DGP1?
    for (r in 1:R){
        if (DGP == "DGP2")
        { d0 <- DGP2(N, K, rho, K0, K.red, b1, b2, p1, p2) }
        else if (DGP == "DGP3") 
        {d0 <- DGP3(N, K, rho, K0, K.red, b1, b2, p1, p2)}
        else if (DGP == "DGP4") 
        {d0 <- DGP4(N, K, rho, K0, K.red, b1, b2, p1, p2, dist(K))}
        
        
        b <- est(d0$y, d0$X, d0$Z, d0$b0) 
        
        B[, , r] <- b
        if (r%%100 == 1){     print(r);     print( colSums( (b - d0$b0)^2 ) ) }
    }
    
    summ <- summ.stat(B, d0$b0)
    print(summ)
    return(list(summ = summ, B = B) )
}


bhat.oracle <- function(y, X, Z, b0){
  # the oracle estimator of beta
  # when the location is known
  K <- length(b0)
  bhat <- rep(0, K)
  oracle.set <- ( abs(b0) > 0.1 )
  X.oracle <- X[, oracle.set]
  PZ <- Z %*% solve(crossprod(Z)) %*% t(Z)
  bhat.o <- solve( t(X.oracle) %*% PZ  %*% X.oracle ) %*% ( t(X.oracle) %*% PZ  %*% y )
  bhat[oracle.set] <- bhat.o
  return(bhat)
}

#################################

est <- function(y, X, Z, b0)
  # combined estimation
  # it includes the four estimators,
  # GMM-Lasso-AIC, BIC, standard TSLS and oracle
{
  bhat.lars <- gmm_s(y, X, Z)
  bhat.las  <- bhat.lars$bhat
  
  bhat.tsls <- LIML(y, X, Z, "tsls")
  # bhat.liml <- LIML(y, X, Z, "liml") # liml was commented out
  
  #### BIC ######
  mB <- mIC(y, X, Z, bhat.las, "BIC")
  choice.mB <- which.min(mB)
  mA <- mIC(y, X, Z, bhat.las, "AIC")
  choice.mA <- which.min(mA)
  ###############
  
  bhat.o <- bhat.oracle(y, X, Z, b0)
  
  bhat <- cbind( bhat.las[, c(choice.mA, choice.mB )], bhat.lars$bhat.cv, bhat.tsls, bhat.o)
  return(bhat)
}



##############################

two.equation.est <- function(k0, y, X, Z){
    # this function is only useful for inference.
    # it is irrelevant for the estimation, the current paper on ER
    
    bhat.s <- gmm_s(y, X, Z) # structural equation
    bhat.r <- gmm_s(X[,k0], X[, -k0], Z) # reduced-form equation
    
    
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
    # projection matrix
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
