rm(list = ls())
setwd("H://Dropbox/code/aug_OLS_R")
source("DGP.R")
source("LIML.R")
source("gmm.R")
library(lars)
library(svd)
set.seed(77)
#######################

summ.stat <- function(B, b0){
  biasS <- colSums( (apply( B, c(1,2), mean) - b0)^2 )
  vars <- colSums( apply(B, c(1,2), var) )
  MSE  <- biasS + vars
  tab <- rbind( biasS, vars, MSE)
  colnames(tab) <- c("AIC", "BIC", "Oracle", "GMM", "LIML")
  return(tab )
}
############################s
#############################
N <- 300
K <- 40
K0 <- 10
K.red <- 1
rho <- 0.3

b1 <- 1
b2 <- 0.0

R <- 20
#############
B <- array(0, dim = c(K, 5, R) )

#############

for (r in 1:R){
  d0 <- DGP1(N, K, rho, K0, K.red, b1, b2)
  
  b <- est(d0$y, d0$X, d0$Z, d0$b0)
  print(r)
  B[, , r] <- b
  print( colSums( (b - d0$b0)^2 ) )
}

summ <- summ.stat(B, d0$b0)
print(summ)


###################

