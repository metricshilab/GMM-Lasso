# 

K <- 80
p1 <- 1
AA <- matrix( 0, ncol = K, nrow = K)
mult <- sqrt(1/3)
diag(AA) <- mult * p1
AA[1, K] <- mult * p1
AA[K, 1] <- mult * p1
AA[seq(from = 2, to = K*K, by = K+1)] <- mult * p1
AA[seq(from = K+1, to = K*K, by = K+1)] <- mult * p1

eig <- eigen(crossprod(AA), only.values = T)$values
print(min(eig))