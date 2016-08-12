
############################
MSE <- function(bbb, b.true) {
  bias <- colMeans( bbb[, 1:2] - b.true )
  var  <- apply( bbb[, 1:2], 2, var )
  mse  <- bias^2  + var
  return(data.frame( bias = bias, var = var, mse = mse) )             
}

##################
workhorse <- function(N,K, K0, K.red, rho, b1, b2, p1, k0, DGP, R = 20){
  b.true <- b1
  p2 <- 1 - rho - p1

  B <- matrix(NA, ncol = 6, nrow = R)
  test <- rep(0, R)
  
  r <- 1
  while (r <= R){
  if (DGP == "DGP2") {
    d0 <- DGP2(N, K, rho, K0, K.red, b1, b2, p1, p2)  
  } else if (DGP == "DGP4") {
    d0 <- DGP4(N, K, rho, K0, K.red, b1, b2, p1, p2)
    }

    bb <- two.equation.est_iter(k0, d0$y, d0$X, d0$Z)
    
    bbb <- stat(k0, d0$y, d0$X, d0$Z, bb, display = F)
    
    B[r, ] <- bbb
    test[r] <- (B[r,3] < b.true) & (b.true < B[r,4] )
    ratio <- sum(test)/r
    if (r%%100 == 1) { cat("r = ", r, "ratio = ", ratio, "\n" ) }
    
    r <- r + 1
    if( bb$flag == 0 ) { r <- r-1 }
  }
  colnames(B) <- c("b.init", "b.corr", "CI.L", "CI.U", "sig.e", "tau")
  
  beta <- B[, c(1:2, 5:6)]
  coverage <- mean(test)  
  summ.stat <- MSE(B, b.true) 
  return( list(  b = beta, coverage = coverage, summ = summ.stat) )
}
##########33
