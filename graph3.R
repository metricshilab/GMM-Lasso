rm(list = ls())
library(ggplot2)
library(reshape2)
library(plyr)

load("DGP2 KK 20, 40, 80 b1 1 b2 0 p1 0.3 rho 0 .Rdata")

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


pc2 <- pc1+ ylab("")+labs(title = paste("b1 =", b1, ", b2 =", b2, ", gamma2 = ", rho))

pc2 <- pc2 + theme_bw() + theme(panel.grid.minor = element_blank())
pc2 <- pc2 + scale_fill_manual(values=c("black", "grey"))

print(pc2)


