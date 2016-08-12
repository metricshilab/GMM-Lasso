# setwd("/estimation_plot")

library(ggplot2)
library(reshape2)
library(plyr)
setwd("results_est/")
####
filenames <- list.files()

myplot <- function(file){
  load(file)
  
  # pc2 <- plotgraph(result, b1, b2, p1, rho)  
  p <- ggplot(data = result, aes(x = variable, fill = measure, weight = value))
  
  pc1 <- p + geom_bar(position = "stack")
  pc2 <- pc1 + facet_grid(N~K, scales = "free", labeller = label_both)
  pc2 <- pc2 + xlab("Estimators") + ylab("squared bias and variance") +
    labs(title = paste("b1 =", b1, ", b2 =", b2, ", gamma1 =", p1, ", 
                           gamma2 = ", rho)             )
  # pc3 <- pc2 + theme_grey()
  # print(pc2)  
  ggsave(pc2, file = paste(filename,".png" ), 
         width = 8.22,  height = 6.54, unit = "in")
}

l_ply( .data = filenames, .fun = myplot )

