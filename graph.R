#################################################
plotgraph.bias <- function(result, b1, b2, p1, rho){

  matrixA <-  result[[1]]  #; browser()

  for (i in 2:length(result) ){
    matrixA <- rbind(matrixA, result[[i]] )
  }


  KKK <- data.frame(N      = as.factor(matrixA[,1]) ,
                    K      = as.factor(matrixA[,2]) ,
                    b.init = matrixA[,3], 
                    b.corr = matrixA[,4]
                    )

  
  KKKK <- melt(KKK, id.vars = c("N", "K"), 
              measure.vars= c("b.init", "b.corr") )
######################
  p <- ggplot(data = KKKK, aes(x = value, fill = variable, alpha = 0.05) )
  pc1 <- p + geom_density() + geom_vline( xintercept = b1)
  pc2 <- pc1 + facet_grid(N~K, scales = "free", labeller = label_both)
  pc2 <- pc2 + xlab("")  + labs(title = paste("b1 =", b1, ", b2 =", b2) )
  print(pc2)
  
  return(pc2)            
}
#####################################################

plotgraph.normal <- function(result, b1, b2, p1, rho){

  matrixA <-  result[[1]]  #; browser()
  
  for (i in 2:length(result) ){
    matrixA <- rbind(matrixA, result[[i]] )
  }
  

  KKK <- data.frame(N      = as.factor(matrixA[,1]) ,
                    K      = as.factor(matrixA[,2]) ,
                    # standardize point estimate
                    b.corr.std = ( matrixA[,4] - b1 ) * matrixA[ ,5]/ matrixA[,6] 
  )
  
  
  KKKK <- melt(KKK, id.vars = c("N", "K"), 
               measure.vars= c("b.corr.std") )
  base <- seq(-3, 3, length.out = 200 )
  data.normal <- data.frame( base = base,  y.base = dnorm(base) )
  ######################
  p <- ggplot(data = KKKK, aes(x = value) )
  pc0 <- p +geom_density(fill = "red", alpha = 0.5, color = "red") 
  pc1 <- pc0  + geom_line( aes(x = base, y = y.base ), data = data.normal, size = 1)  
  pc2 <- pc1 + facet_grid(N~K, scales = "free", labeller = label_both)
  pc2 <- pc2 + xlab("")  + labs(title = paste("b1 =", b1, ", b2 =", b2) )  + xlim(c(-3,3))
  print(pc2)
  
  return(pc2)            
}
