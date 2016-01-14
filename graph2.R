rm(list = ls())
library(ggplot2)
library(reshape2)
setwd("H://Dropbox/code/aug_OLS_R")
load("H:/Dropbox/code/aug_OLS_R/K.red 1 b2 0 0.3 0.3 .Rdata")

convert <- function(A){
  # A is an object in the list
  A <- data.frame(A)
  AA <- melt(A, id = 
}

resultA <- list()
resultA[1:3] <- result[c(1,2,2)]
resultA[4:9] <- result[3:8]
resultA[[3]][,] <- 0

matrixA <- t( resultA[[1]] )
for (i in 2:9){
  matrixA <- rbind(matrixA, t(resultA[[i]]) )
}
AA[, 4:6] <- matrixA

AAA <- AA[1:5, 1:5]
AAAA<- melt(AAA, id = c("N", "K", "Estimator"), variable.name = "measure")

p <- ggplot(data = AAAA, aes(x = Estimator, fill = measure, weight = value))
p1 <- p + geom_bar(position = "stack")
print(p1)
