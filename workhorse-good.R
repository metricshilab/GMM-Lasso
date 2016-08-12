rm(list = ls())
setwd("H://Dropbox/code/aug_OLS_R")
source("DGP.R")
source("gmm_pen.R")
source("gmm_iter_pen.R")
source("funcs_inference.R")

library(lars)
library(svd)
library(ggplot2)
library(reshape2)
library(penalized)

# set.seed(99)

############################
MSE <- function(bbb, b.true) {
  bias <- colMeans( bbb[, 1:2] - b.true )
  var  <- apply( bbb[, 1:2], 2, var )
  mse  <- bias^2  + var
  return(data.frame( bias = bias, var = var, mse = mse) )             
}

##################

N <- 200
K <- 20
K0 <- 10
K.red <- 1.1
rho <- 0.3

b1 <- 1
b2 <- 0.0


p1 <- 0.3
p2 <- 1 - rho - p1
k0 <- 1
b.true <- b1


R <- 20
B <- matrix(NA, ncol = 5, nrow = R)
test <- rep(0, R)

r <- 1
while (r <= R){

  d0 <- DGP4(N, K, rho, K0, K.red, b1, b2, p1, p2)  
  # bb <- two.equation.est(k0, d0$y, d0$X, d0$Z)
  bb <- two.equation.est_iter(k0, d0$y, d0$X, d0$Z)
  
  bbb <- stat(k0, d0$y, d0$X, d0$Z, bb, display = F)
  
  B[r, ] <- bbb
  test[r] <- (B[r,3] < b.true) & (b.true < B[r,4] )
  ratio <- sum(test)/r
  cat("r = ", r, "ratio = ", ratio, "\n" )
  
  r <- r + 1
  if( bb$flag == 0 ) { r <- r-1 }
}
colnames(B) <- c("b.init", "b.corr", "CI.L", "CI.U", "sig.e")

print(B)
print(mean(test))


print( MSE(B, b.true) )
