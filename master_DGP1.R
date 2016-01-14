############################
summ.stat <- function(B, b0){
    # it is a function called in "master1"
    # calculate the summary statistics for MSE, squared-bias, and variance
    
  biasS <- colSums( (apply( B, c(1,2), mean) - b0)^2 )
  vars <- colSums( apply(B, c(1,2), var) )
  MSE  <- biasS + vars
  tab <-  rbind( biasS, vars, MSE) 
  colnames(tab) <- c("AIC", "BIC", "GMM", "Infeasible")
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

  #############
  B <- array(0, dim = c(K, 4, R) )

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

