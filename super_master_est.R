# super master file to run the simulations for gmm lasso
# also plot the graphs

# last update 2014-11-14

rm(list = ls())

source("DGP.R")
source("estimation.R")
source("graph.R")
source("cv.lars.2sls.R")

package.list <- c("lars", "svd", "ggplot2", "reshape2", "plyr")
# install.packages(package.list)
library(lars)
library(svd)
library(ggplot2)
library(reshape2)
library(plyr)
library(doParallel)
library(foreach)

ptm <- proc.time()

set.seed(88)

supermaster <- function(b1, b2, rho, DGP){
  # specify the parameters   
  R <- 1000 # number of replications
  p1 <- 0.3
  
  NN <- c(200, 400) # sample size
  KK <- c( 20, 40, 80)    # the total number of regressors
  K0 <- 10               # the number of "large" coefficients
  K.red <- 1.1           # add some redudant IVs to for over-identification
  
  p2 <- 1 - rho - p1
  filename <- paste( "R", R, DGP,"KK", toString(KK), "b1", b1, "b2", b2, "p1", p1, "rho", rho, ".Rdata", sep = "")
  
  #### start the loop
  r <- 0
  
  result <- list()
  resultB <- list()
  
  for (N in NN){
    for ( K in KK ){
      cat("K = ", K, "N = ", N, "\n")
      r <- r + 1
      # if ( K * K.red * 1.5 < N ){
      # "master1" is the key function that does the actual job
      # "master1" is a function in "estimation.R"
      res <- master1(N,K, K0, K.red, rho, b1, b2, p1, p2, R, DGP)
      result[[r]] <- data.frame(N = N, K = K, MSE = res$summ)
      resultB[[r]] <- res$B
      print(r) 
    }
  }  
  
  result <- ldply(result)
  names(result) <- c("N", "K", "AIC", "BIC", "CV", "2SLS", "IFB")
  result$measure <- c("biasSQ", "var", "MSE")
  result <- melt(result, id = c("N", "K", "measure") )
  result <- result[result$measure != "MSE" & result$K != 20, ]
  
  print(filename)
  save(N,K, K0, K.red, rho, b1, b2, p1, R, DGP, result, filename, file = filename)
  
  ### plot
#   pc2 <- plotgraph(result, b1, b2, p1, rho)  
#   print(pc2)  
#   ggsave(pc2, file = paste(filename,".eps" ), 
#          width = 8.22,  height = 6.54, unit = "in")
}
####
# in the simulation I only used two DGPs, DGP2 and DGP4, with different parameters.

parameter.meta <- data.frame( b1 = c(rep(1,4), rep(-1,4), c(1,-1) ) ,
                              b2 = c(0, 0, 0.01, 0.01, 0, 0, -0.01, -0.01, 0, 0), 
                              rho = c( rep(0.3, 8), rep(0, 2) ),
                              dgpID = c( rep(c("DGP2", "DGP4"), 4), "DGP2", "DGP2" ) )

registerDoParallel()

l_ply(.data = 1:10, 
      .fun = function(i) with(parameter.meta, supermaster(b1[i], b2[i], rho[i], dgpID[i]) ) ,
      .parallel = TRUE, 
      .paropts = list( .packages = package.list, .export = ls(envir=globalenv() ) )
)
# strange enough. The order of the arguments in ".paropts" seem to matter
# I will check this issue later

print(proc.time() - ptm )
# supermaster(1, 0, "DGP2")
# supermaster(1, 0, "DGP4")

# supermaster(1, 0.01, "DGP2")
# supermaster(1, 0.01, "DGP4")

# supermaster(-1, 0, "DGP2")
# supermaster(-1, 0, "DGP4")

# supermaster(-1, -0.01, "DGP2")
# supermaster(-1, -0.01, "DGP4")

# two complementary graphs
# supermaster(1,0, "DGP2")
# supermaster(-1,0, "DGP2")