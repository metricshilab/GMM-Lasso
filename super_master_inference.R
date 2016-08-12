# super master file to draw the graphs using the data generated from super master 1
# last update 2014-11-12

rm(list = ls())

source("DGP.R")
source("gmm_pen.R")
source("gmm_iter_pen.R")
source("funcs_inference.R")
source("plot_workhorse.R")
source("graph.R")

library(lars)
library(svd)
library(ggplot2)
library(reshape2)
library(penalized)

### workhorse(200,20, 10, 1.1, 0.3, 1, 0, .3, 1, "DGP2") 
#######################

supermaster2 <- function(b1, b2, k0, rho, DGP, R = 500) {

  p1 <- 0.3
  
  result <- list()
  resultB <- list()
  result.coverage <- numeric()
  
  
  NN <- c(200, 400, 800, 1600)
  KK <- c(20, 40, 80)
  K0 <- 10
  K.red <- 1.1


  filename <- paste("inference", DGP,"KK", toString(KK), "b1", b1, "b2", b2, "p1", p1, "rho", rho, ".Rdata")
  #########################
  
  r <- 0
  
  for (N in NN){
    for ( K in KK ){
      cat("K = ", K, "N = ", N, "\n")
      r <- r + 1
      if ( K * K.red * 1.5 < N ){
            
        res <- workhorse(N,K, K0, K.red, rho, b1, b2, p1, k0, DGP, R)
        result[[r]] <- res$summ
        dimension <- matrix( c(N, K), byrow = T, ncol = 2, nrow = R ) 
        resultB[[r]] <- cbind( dimension, res$b )
        result.coverage[r] <- res$coverage
        print(r)
      }
      else {
        result[[r]] <- result[[r-1]]
        resultB[[r]]<- resultB[[r-1]]
        result.coverage[r] <- 0
        result[[r]][,] <- 0
        resultB[[r]][,,] <- 0
        }
    }
  }

  kaka <- list( B = resultB, mse = result, coverage = result.coverage ) 
  
  save(N,K, K0, K.red, rho, b1, b2, p1, R, DGP, kaka, file = filename)

  pc2 <- plotgraph.bias(resultB, b1, b2, p1, rho)  
  print(pc2)  
  ggsave(pc2, file = paste(filename,".png" ), width = 8.22,  height = 6.54, unit = "in")
  
  pc3 <- plotgraph.normal(resultB, b1, b2, p1, rho)
  ggsave(pc3, file = paste("distribution", filename,".png" ), width = 8.22,  height = 6.54, unit = "in")
}
####################################################


# supermaster(1, 0, "DGP2")
# supermaster(1, 0, "DGP4")

# supermaster(1, 0.01, "DGP2")
# supermaster(1, 0.01, "DGP4")

# supermaster(-1, 0, "DGP2")
# supermaster(-1, 0, "DGP4")

# supermaster(-1, -0.01, "DGP2")
# supermaster(-1, -0.01, "DGP4")

# two complementary graphs
a <- supermaster2(b1 = -1, b2 = 0, k0 = 1, rho = .3, "DGP4", R = 500)
# supermaster2(-1,0, "DGP2")

