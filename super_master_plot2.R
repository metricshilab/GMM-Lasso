

library(ggplot2)
library(reshape2)
library(plyr)
#setwd("results_est/")
####
filenames <- list.files()

myplot <- function(file){
  
  load(file)
  
  result <- ldply(result)
  
  ID <- data.frame(N = rep(c(20,40,80)*10, each = 9), 
                   K = rep( rep( c(200, 400, 800)/10, each = 3 ) , 3 ),
                   Estimator= rep(c("biasSQ", "var", "MSE"), 9)
  )
  
  res <- cbind(ID, result)
  
  res1<- melt(res, id = c("N", "K", "Estimator"), 
              measure.vars = c("AIC", "BIC", "GMM", "Infeasible"))
  res2 <- res1[ res1$Estimator != "MSE", ]
  
  names(res2) <- c("N", "K", "measure", "estimator", "value")
  
  
  
  ######### plot 
  p <- ggplot(data = res2, aes(x = estimator, y = value, fill = measure ) )
  pc1 <- p + geom_bar(stat = "identity") + facet_grid( N~K , scales = "free", labeller = label_both)
  
  if (b2 == 0){ 
    sparse.type = "Exactly Sparse"
  } else {
    sparse.type = "Approximately Sparse"
  }
    
  pc2 <- pc1+ ylab("")+labs(title = sparse.type )+
    xlab( paste("parameters: b1 =", b1, ", b2 =", b2, ", gamma2 = ", rho))
  
  pc2 <- pc2 + theme_bw() + theme(panel.grid.minor = element_blank(), strip.background = element_blank())
  pc2 <- pc2 + scale_fill_manual(values=c("black", "grey"))
  

  ggsave(pc2, file = paste(file,".eps" ), 
         width = 8.22,  height = 6.54, unit = "in")
}

l_ply( .data = filenames, .fun = myplot )
