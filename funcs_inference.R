# for each model, calculate the relevant quantitiles

lambda1 <- function(model){
  lambda <- model@lambda1
  }

######################################

gmm_s_iter_reduced_1 <- function(y, X, Z, display = F, tuning = F){
  
  browser() 
  
  flag <- 1
  
  sig.e1 <- 1
  sig.e2 <- 0
  Tol <- 0.05
  m <- 0
  
  PZ <- Proj(Z)
  
  while ( abs(sig.e1 - sig.e2 ) > Tol ){
    m <- m + 1
    
    b1 <- lasso_lambda(y/sig.e1, X/sig.e1, Z, lambda)
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

##################################3
lasso_lambda<- function(y, X, lambda){
#   
#   D <- dim(X)
#   N <- D[1]
#   K <- D[2] 
#   
#   PZ <- solve( t(Z) %*% Z )
#   svd.PZ <- svd(PZ)
#   HProj <- diag(sqrt(svd.PZ$d)) %*% t(svd.PZ$v) %*% t(Z)
#   
#   PX = HProj %*% X
#   Py = HProj %*% y
#   
#   # reg <- lars(PX, Py, type = "lasso", normalize = F, intercept = F)
  
  reg <- penalized( response = Py, penalized =  PX , 
                    unpenalized = ~0, 
                    lambda1 = lambda, 
                    model = "linear")
  beta <- coef(reg)
  res <- residuals(reg)
  return(beta = beta, res = res)                  
  
  }