gmm_s_AIC <- function(y, X, Z, tune = F){
  # _s indicates the structural equation
  # use "lars" to compute the path of the lasso solutions
  
  bhat <- gmm_s(y, X, Z)
  #####################
  mAIC <- mIC(y, X, Z, bhat, "AIC")
  if (tune == F){
    choice <- which.min(mAIC)
  } else {
    choice <-  min( which.min(mAIC) * 3, length(mAIC) )
    print(choice)
    }
  b <- bhat[,  choice]
  return(b) 
}

#####################################


gmm_s_iter <- function(y, X, Z, display = F){
  
  flag <- 1
  
  sig.e1 <- 1
  sig.e2 <- 0
  Tol <- 0.05
  m <- 0
  while ( abs(sig.e1 - sig.e2 ) > Tol ){
    m <- m + 1

    b1 <- gmm_s_AIC(y/sig.e1, X/sig.e1, Z)
    e2 <- as.vector( y - X%*%b1 ) 
    sig.e2 <- sqrt( sd(e2) )# update sig.e1
    
    b2 <- gmm_s_AIC(y/sig.e2, X/sig.e2, Z)
    e1 <- as.vector( y - X%*%b2 )
    sig.e1 <- sqrt( sd(e1) ) # update sig.e2

    if (display == T) { cat( "diff = ",  abs (sig.e1 - sig.e2), "\n" ) }
    if (m > 20) {  flag = 0; break }
  }
  if (display == T) {   cat("sig.e = ", sig.e1, "\n") }
  return(list( b2 = b2, flag = flag ) )
}

#############################


gmm_s_iter_r <- function(y, X, Z, display = F, tuning = F){
  
  flag <- 1
  
  sig.e1 <- 1
  sig.e2 <- 0
  Tol <- 0.05
  m <- 0
  
  PZ <- Proj(Z)
  
  while ( abs(sig.e1 - sig.e2 ) > Tol ){
    m <- m + 1
    
    b1 <- gmm_s_AIC(y/sig.e1, X/sig.e1, Z)
    e2 <- as.vector( y - X%*%b1 ) 
    # browser()
    sig.e2 <- as.vector( sqrt( t(e2)  %*% PZ %*% e2 ) )# update sig.e1
    
    b2 <- gmm_s_AIC(y/sig.e2, X/sig.e2, Z)
    e1 <- as.vector( y - X%*%b2 )
    sig.e1 <- as.vector( sqrt(  t(e1)  %*% PZ %*% e1) ) # update sig.e2
    
    if (display == T) { cat( "diff = ",  abs (sig.e1 - sig.e2), "\n" ) }
    if (m > 20) {  flag <- 0; break }
  }
  if (display == T) {   cat("sig.e = ", sig.e1, "\n") }
  return(list(b2 = b2, flag = flag) )
}

#############################
#############################

two.equation.est_iter <- function(k0, y, X, Z){
  
  bhat.s <- gmm_s_iter(y, X, Z)
  bhat.r <- gmm_s_iter_r(X[,k0], X[, -k0], Z)
  flag   <-  bhat.s$flag * bhat.r$flag
  bhat <- list(s = bhat.s$b2, r = bhat.r$b2, flag = flag)
}
########################################
