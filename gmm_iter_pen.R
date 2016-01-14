# the functions in this script

# gmm_s_AIC: esitmation of the structural equation using AIC
# gmm_s_iter_s: estimation of the structural equation using the iteration
#               mainly for a better estimation of the variance
# gmm_pen_reduce: estimation of the reduce-form using "penalized".
#                 the iteration procedure is my creation
# two.equaiton.est_iter: estimation and return the betas in both the structural 
#                        equation and the reduced-form equation.
# b.compare: the function of choice along the path of lambdas, which generate a path of betas


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
gmm_s_iter_s <- function(y, X, Z, display = F){
  # iterative estimation of the structural model
  
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

    # if (display == T) { cat( "diff = ",  abs (sig.e1 - sig.e2), "\n" ) }
    if (m > 20) {  flag = 0; break }
  }
  if (display == T) {   cat("sig.e = ", sig.e1, "\n") }
  return(list( b2 = b2, flag = flag ) )
}
#############################
gmm_pen_reduce <- function(y, X, Z,lambda, display = F){
  # estimation of the reduced-form using "penalized"
  
  # lambda is a scalar
  # the iteration is not that of Zhang and Zhang. It is mine.
  n <- length(y)
  
  flag <- 1
  
  sig.e1 <- 1
  sig.e2 <- 0
  Tol <- 0.05
  m <- 0
  
  PZ <- Proj(Z)
  
  while ( abs(sig.e1 - sig.e2 ) > Tol ){
    m <- m + 1
    
    b11 <- gmm_s_pen(y/sig.e1, X/sig.e1, Z, lambda)
    b1  <- coef(b11, "all")
    e2 <-   y - X%*%b1 
    sig.e2 <- as.vector( sqrt( t(e2)  %*% PZ %*% e2/n ) )# update sig.e1
    
    b22 <- gmm_s_pen(y/sig.e2, X/sig.e2, Z, lambda)
    b2  <- coef(b22, "all")
    e1 <-   y - X%*%b2 
    sig.e1 <- as.vector( sqrt(  t(e1)  %*% PZ %*% e1 /n) ) # update sig.e2
    
    if (display == T) { cat( "diff = ",  abs (sig.e1 - sig.e2), "\n" ) }
    if (m > 20) {  flag <- 0; break }
  }
  # if (display == T) {   cat("sig.e = ", sig.e1, "\n") }
  return(list(b2 = b2, flag = flag) )
}
#############################
two.equation.est_iter <- function(k0, y, X, Z){
  K <- ncol(X) - 1

  
  bhat.s <- gmm_s_iter_s(y, X, Z)
  
  
  lambda <- get.lambda(X[,k0], X[,-k0], Z, round(K/2))
  
  
  # lambda <- seq(from = .001, to = 500, length.out = 2)
  lambda.length <- length(lambda)
  
  bhat.r.b <- matrix(0, ncol = lambda.length, nrow = K )
  bhat.r.f <- rep(0, lambda.length)
  # bhat.mIC <- rep(0, lambda.length)
  
  
  for ( mm in 1:lambda.length ){
    lambda.t <- lambda[mm]
    bhat.rr <- gmm_pen_reduce(X[,k0], X[, -k0], Z, lambda.t )
    
    bhat.r.b[, mm ]<- bhat.rr$b2 # coef(bhat.rr, "all") 
    bhat.r.f[mm] <- bhat.rr$flag
  }
  bmIC <- mIC(X[,k0], X[, -k0], Z, bhat.r.b, "AIC")

  flag   <-  1 # bhat.s$flag * bhat.r$flag

  bhat.r <- b.compare( X[,k0], X[, -k0], Z, bhat.r.b, lambda, which.min(bmIC), 4  )
  bhat <- list(s = bhat.s$b2, r = bhat.r, flag = flag)
}
########################################
b.compare <- function(y, X, Z, b, lambda, choice, multiplier){
  # multiplier: a constant > 1 to allow some relaxation of the variance
  #             used to trade some variance for a better bias.
  e <- y - X %*% b
  square.root.uPu <- sqrt( diag( t( e ) %*% Z %*% solve( t(Z) %*% Z ) %*% t(Z) %*% e ) )
  uPX <- abs( diag( t( e ) %*% Z %*% solve( t(Z) %*% Z ) %*% t(Z) %*% X ) ) 

  var  <- square.root.uPu/uPX
  bias <- lambda/ square.root.uPu


  index <- ( var <= var[choice] * multiplier )
  index.b <- which.min( bias[index] )
  final.b <- b[, index.b]

  return(final.b)
  }