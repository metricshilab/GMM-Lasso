# cv.lars.2sls

cv.lars.2sls <- function (x, y, z, K, index, trace = FALSE, plot.it = TRUE, 
                          se = TRUE, type = c("lasso", "lar", "forward.stagewise", 
                          "stepwise"), mode = c("fraction"), ...) 
{
  # x, y, z are the original data
  # will do the transformation inside the function
  
  type = match.arg(type)
  if (missing(mode)) {
    mode = switch(type, lasso = "fraction", lar = "step", 
                  forward.stagewise = "fraction", stepwise = "step")
  }
  else mode = match.arg(mode)
  
  all.folds <- cv.folds(length(y), K)
  
  if (missing(index)) {
    index = seq(from = 0, to = 1, length = 100)
  }
  
  residmat <- matrix(0, length(index), K)
  for (i in seq(K)) {
    omit <- all.folds[[i]]
    
    x.omit <- x[-omit, ]
    y.omit <- y[-omit]
    z.omit <- z[-omit, ]
    
    PXy.omit <- tsls.transform(y = y.omit, X = x.omit, Z = z.omit)
    PX.omit  <- PXy.omit$PX 
    Py.omit  <- PXy.omit$Py 
    
    reg <- lars( x = PX.omit, y = Py.omit, trace = trace, type = type, ...)
    
    # need to fix the bug about the criterion of fit. 
    # also need to check the old code of AIC and BIC.
    
    fit <- predict(reg, newx = x.omit, mode = mode, s = index)$fit
    
    if (length(omit) == 1) 
      fit <- matrix(fit, nrow = 1)
    
    Pfit <- tsls.transform( y = y.omit, X = fit, Z = z.omit)$PX # transform "fit" into the desirable shape

    residmat[, i] <- apply(( as.vector(Py.omit) - Pfit)^2, 2, mean)
    if (trace) cat("\n CV Fold", i, "\n\n")
  }
  cv <- apply(residmat, 1, mean)
  cv.error <- sqrt(apply(residmat, 1, var)/K)
  object <- list(index = index, cv = cv, cv.error = cv.error, 
                 mode = mode)
  if (plot.it) 
    plotCVLars(object, se = se)
  
  return(object)
}
