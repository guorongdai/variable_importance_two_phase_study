suppressMessages({
  
  if (!require(randomForest)) install.packages('randomForest', repos =
                                                 'https://cran.revolutionanalytics.com/')
  
  #  if (!require(quantregForest)) install.packages('quantregForest', repos =
  #                                                   'https://cran.revolutionanalytics.com/')
  
  #  if (!require(quantreg)) install.packages('quantreg', repos =
  #                                             'https://cran.revolutionanalytics.com/')
  
  if (!require(gbm)) install.packages('gbm', repos =
                                        'https://cran.revolutionanalytics.com/')
  
  if (!require(mgcv)) install.packages('mgcv', repos =
                                         'https://cran.revolutionanalytics.com/')
  
  if (!require(xgboost)) install.packages('xgboost', repos =
                                            'https://cran.revolutionanalytics.com/')
  
  if (!require(np)) install.packages('np', repos =
                                       'https://cran.revolutionanalytics.com/')
  
  if (!require(qgam)) install.packages('qgam', repos =
                                         'https://cran.revolutionanalytics.com/')                                    
  
  library(randomForest)
  #  library(quantregForest)
  #  library(quantreg)
  library(gbm)
  library(mgcv)
  library(xgboost)
  library(np)
  library(qgam)
  
})

MVIM = function(Y, X, Z, K, method1, tau = NULL, loss, depth1 = NULL)
{
  
  if(loss == "quadratic") M = 1
  if(loss == "check") M = length(tau)
  
  U = cbind(X, Z)
  
  ghat = CM.fit(U, Y, K, M, method1, tau, depth1)
  hhat = CM.fit(X, Y, K, M, method1, tau, depth1)
  
  if(loss == "quadratic")
  {
    
    vehat = as.matrix( quadratic(Y - ghat[, 1]) )
    etahat = as.matrix( quadratic(Y - hhat[, 1]) )
    
  }
  
  if(loss == "check")
  {
    
    vehat = multiple.check(Y, ghat, tau)
    etahat = multiple.check(Y, hhat, tau)
    
  }
  
  thetahatm1 = colMeans(vehat)
  thetahatm2 = colMeans(etahat) 
  
  thetahatm = thetahatm1 / thetahatm2
  thetahat = mean(thetahatm)
  
  result = list( "thetahatm" = thetahatm, "thetahat" = thetahat)
  
  return(result)
  
}

#########################################################################################################
# Functions for Poisson sampling
MVIM.TPD.B = function(Y, X, Z, R, pi, K, method1, method2, tau = NULL, loss, depth1 = NULL, depth2 = NULL)
{
  
  if(loss == "quadratic") M = 1
  if(loss == "check") M = length(tau)
  
  N = length(Y)
  
  W = cbind(as.matrix(Y), X)
  U = cbind(X, Z)
  S = R / pi
  
  gtilde = CM.fit.B(U, Y, R, pi, K, M, method1, tau, depth1)
  hhat = CM.fit(X, Y, K, M, method1, tau, depth1)
  
  if(loss == "quadratic")
  {
    
    vetilde = as.matrix( quadratic(Y - gtilde[, 1]) )
    etahat = as.matrix( quadratic(Y - hhat[, 1]) )
    
  }
  
  if(loss == "check")
  {
    
    vetilde = multiple.check(Y, gtilde, tau)
    etahat = multiple.check(Y, hhat, tau)
    
  }
  
  mutilde = NF.Est.B(Y, X, Z, R, pi, K, M, method1, method2, loss, tau, depth1, depth2)
  
  ############################################################
  # our estimator and its variance 
  thetatildem1 = colMeans(mutilde + (vetilde - mutilde) * S )
  thetahatm2 = colMeans(etahat) 
  thetatildem = thetatildem1 / thetahatm2
  thetatilde = mean(thetatildem)
  
  phihatm = t( t(mutilde + (vetilde - mutilde) * S) / thetahatm2 - thetatildem1 * t(etahat) / thetahatm2 ^ 2 )
  phihat = rowMeans(phihatm)
  lambdatildem = sqrt( colVar(phihatm) )
  lambdatilde = sqrt( var(phihat) ) 
  ############################################################
  
  ############################################################
  # zero-imputation estimator and its variance 
  thetatildemzero1 = colMeans(vetilde * S)
  thetahatmzero2 = colMeans(etahat)
  thetatildemzero = thetatildemzero1 / thetahatmzero2
  thetatildezero = mean(thetatildemzero)
  
  phihatmzero = t( t(vetilde * S) / thetahatmzero2 - thetatildemzero1 * t(etahat) / thetahatmzero2 ^ 2 )
  phihatzero = rowMeans(phihatmzero)
  lambdatildemzero = sqrt( colVar(phihatmzero) )
  lambdatildezero = sqrt( var(phihatzero) ) 
  ############################################################
  
  ############################################################
  # phase-II-only estimator and its variance
  thetatildeOm1 = colMeans(vetilde * S)
  thetatildeOm2 = colMeans(etahat * S)
  thetatildeOm = thetatildeOm1 / thetatildeOm2
  thetatildeO = mean(thetatildeOm)
  
  phihatOm = t( t(vetilde * S) / thetatildeOm2 - thetatildeOm1 * t(etahat * S) / thetatildeOm2 ^ 2 )
  phihatO = rowMeans(phihatOm)
  lambdatildeOm = sqrt( colVar(phihatOm) )
  lambdatildeO = sqrt( var(phihatO) )
  ############################################################
  
  result = list( "thetatildem" = thetatildem, "sdtildem" = lambdatildem / sqrt(N),
                 "thetatilde" = thetatilde,  "sdtilde" = lambdatilde / sqrt(N),
                 "thetatildemzero" = thetatildemzero, "sdtildemzero" = lambdatildemzero / sqrt(N),
                 "thetatildezero" = thetatildezero, "sdtildezero" = lambdatildezero / sqrt(N),
                 "thetatildeOm" = thetatildeOm, "sdtildeOm" = lambdatildeOm / sqrt(N),
                 "thetatildeO" = thetatildeO, "sdtildeO" = lambdatildeO / sqrt(N) )
  
  
  return(result)
  
}


CM.fit.B = function(X, Y, R, pi, K, M, method, tau = NULL, depth = NULL)
{
  
  n = length(Y)
  nK = floor(n / K)
  
  ob = which(R == 1)
  
  result = matrix(0, n, M)
  
  for(k in 1 : (K - 1))
  {
    
    ind = (k - 1) * nK + (1 : nK)
    ind.pre = intersect(ind, ob)
    ind.fit = intersect( (1 : n)[-ind], ob )
    
    Y.fit = Y[ind.fit]
    X.fit = as.matrix(X[ind.fit, ])
    pi.fit = pi[ind.fit]
    X.predict = as.matrix(X[ind.pre, ])
    
    result[ind.pre, ] = CM.Est.B(Y.fit, X.fit, 1 / pi.fit, X.predict, method, tau, depth)
    
  }
  
  ind = ( ((K - 1) * nK) + 1 ) : n
  ind.pre = intersect(ind, ob)
  ind.fit = intersect( (1 : n)[-ind], ob )
  
  Y.fit = Y[ind.fit]
  X.fit = as.matrix(X[ind.fit, ])
  pi.fit = pi[ind.fit]
  X.predict = as.matrix(X[ind.pre, ])
  
  result[ind.pre, ] = CM.Est.B(Y.fit, X.fit, 1 / pi.fit, X.predict, method, tau, depth)
  
  return(result)
  
}

# Nuisance function estimation with data-splitting
NF.Est.B = function(Y, X, Z, R, pi, K, M, method1, method2, loss, tau, depth1, depth2)
{
  
  n = length(Y)
  nK = floor(n / K)
  
  result = matrix(0, n, M)
  
  W = cbind(as.matrix(Y), X)
  
  for(k in 1 : (K - 1))
  {
    
    ind = (k - 1) * nK + (1 : nK)
    
    Y.fit = Y[-ind]
    X.fit = as.matrix(X[-ind, ])
    Z.fit = as.matrix(Z[-ind, ])
    R.fit = R[-ind]
    pi.fit = pi[-ind]
    W.predict = W[ind, ]
    
    NF.fit = TS.NP.DS.B(Y.fit, X.fit, Z.fit, R.fit, pi.fit, W.predict, M, method1, method2, loss, tau, depth1, depth2)
    
    result[ind, ] = NF.fit
    
  }
  
  ind = ( ((K - 1) * nK) + 1 ) : n
  
  Y.fit = Y[-ind]
  X.fit = as.matrix(X[-ind, ])
  Z.fit = as.matrix(Z[-ind, ])
  R.fit = R[-ind]
  pi.fit = pi[-ind]
  W.predict = W[ind, ]
  
  NF.fit = TS.NP.DS.B(Y.fit, X.fit, Z.fit, R.fit, pi.fit, W.predict, M, method1, method2, loss, tau, depth1, depth2)
  
  result[ind, ] = NF.fit
  
  return(result)
  
}

# Two-stage nonparametric estimator with data-splitting
# method1 is the method used for the first-stage nonparametric estimation
# method2 is the method used for the second-stage nonparametric estimation
TS.NP.DS.B = function(Y, X, Z, R, pi, W.predict, M, method1, method2, loss, tau = NULL, depth1 = NULL, depth2 = NULL)
{
  
  n = length(Y)
  
  n.half = floor(n / 2)
  
  ob = which(R == 1)
  
  W = cbind(as.matrix(Y), X)
  U = cbind(X, Z)
  
  result = matrix(0, nrow(W.predict), M)
  
  for(j in 1 : 2)
  {
    
    if(j == 1)
    {
      ind1 = intersect( 1 : n.half, ob)
      ind2 = intersect( (n.half + 1) : n, ob )
    } else {
      ind2 = intersect( 1 : n.half, ob)
      ind1 = intersect( (n.half + 1) : n, ob )
    }
    
    gtilde.ds = CM.Est.B(Y[ind1], U[ind1, ], 1 / pi[ind1], U[ind2, ], method1, tau, depth1)
    
    if(loss == "quadratic")
    {
      vetilde.ds = quadratic(Y[ind2] - gtilde.ds[, 1])
      result = result + CM.Est.B(vetilde.ds, W[ind2, ], weights = NULL, W.predict, method2, tau = NULL, depth2) / 2
    }
    
    if(loss == "check")
    {
      
      vetilde.ds = multiple.check(Y[ind2], gtilde.ds, tau)
      result = result + apply( vetilde.ds, 2, function(y) CM.Est.B(y, W[ind2, ], weights = NULL, W.predict, method2, tau = NULL, depth2) ) / 2
      
    }
    
  }
  
  return(result)
  
}

CM.Est.B = function(Y.fit, X.fit, weights = NULL, X.predict, method, tau = NULL, depth = NULL)
{
  
  if(!is.null(weights)) weights = weights / sum(weights)
  
  if(method == "GBM.mean")
  {
    
    data.fit = data.frame(Y.fit, X.fit)
    sink("aux")
    model = gbm(Y.fit ~ ., data = data.fit, var.monotone = numeric(ncol(X.fit)),
                distribution = "gaussian", n.trees = 1000, shrinkage = 0.01,
                interaction.depth = depth, bag.fraction = 0.5, train.fraction = 1,
                n.minobsinnode = 10, cv.folds = 5, keep.data = TRUE,
                verbose = F, n.cores = 1)
    best.iter = gbm.perf(model, plot.it = F, method = "cv")
    data.predict = data.frame( numeric(nrow(X.predict)), X.predict )
    result = predict(model, newdata = data.predict, n.trees = best.iter, type = "response")
    sink(NULL)
    
  }
  
  if(method == "GAM.mean")
  {
    
    p = ncol(X.fit)
    
    data.fit = data.frame(Y.fit, X.fit)
    names(data.fit) = c( "Y", paste0("X", 1 : p) )
    X.names = names(data.fit)[-1]
    gam.form = as.formula( paste0("Y ~ ", paste0("s(", X.names, ", bs = 'tp')", collapse = "+")) )
    model = mgcv::gam(gam.form, family = gaussian, data = data.fit, weights = weights,
                      method = "REML", optimizer = c("outer", "bfgs"))
    
    X.predict = data.frame(X.predict)
    names(X.predict) = X.names
    result = predict(model, X.predict, "response")
    
  }
  
  if(method == "GAM.quantile")
  {
    
    p = ncol(X.fit)
    
    data.fit = data.frame(Y.fit, X.fit)
    names(data.fit) = c( "Y", paste0("X", 1 : p) )
    X.names = names(data.fit)[-1]
    gam.form = as.formula( paste0("Y ~ ", paste0("s(", X.names, ", bs = 'tp')", collapse = "+")) )
    sink("aux")
    model = qgam::mqgam( gam.form, data = data.fit, qu = tau,
                         argGam = list(weights = weights, method = "REML", optimizer = c("outer", "bfgs")) )
    sink(NULL)
    X.predict = data.frame(X.predict)
    names(X.predict) = X.names
    Y.predict = qdo(model, tau, predict, newdata = X.predict)
    result = matrix(unlist(Y.predict), ncol = length(tau))
    
  }
  
  return(as.matrix(result))
  
}
#########################################################################################################

#########################################################################################################
# Functions for simple random sampling
# Multiple variable importance measures based two-phase data
# method1 is for the conditional model as well as the first stage of the nuisance estimation
# method2 is for the second stage of the nuisance estimation
# ds indicates whether data-splitting should be used (T) or not (F) in constructing muhat.
# cross = T means when estimating mu(\cdot), we take the average of two data-splitting estimators.
# cross = F means when estimating mu(\cdot), we just take a single data-splitting estimator.
MVIM.TPD = function(Y, X, Z, K, method1, method2, tau = NULL, loss, depth1 = NULL, depth2 = NULL, ds = T, cross = T)
{
  
  if(loss == "quadratic") M = 1
  if(loss == "check") M = length(tau)
  
  N = length(Y)
  n = nrow(Z)
  deltan = n / N
  
  W = cbind(as.matrix(Y), X)
  U = cbind(X[1 : n, ], Z)
  Y.II = Y[1 : n]
  
  ghat = CM.fit(U, Y.II, K, M, method1, tau, depth1)
  hhat = CM.fit(X, Y, K, M, method1, tau, depth1)
  
  if(loss == "quadratic")
  {
    
    vehat = as.matrix( quadratic(Y.II - ghat[, 1]) )
    etahat = as.matrix( quadratic(Y - hhat[, 1]) )
    
  }
  
  if(loss == "check")
  {
    
    vehat = multiple.check(Y.II, ghat, tau)
    etahat = multiple.check(Y, hhat, tau)
    
  }
  
  etahat.II = as.matrix(etahat[1 : n, ])
  
  if(ds == T)
  {
    muhat = NF.Est(Y.II, as.matrix(X[1 : n, ]), Z, W[-(1 : n), ], K, M, method1, method2, loss, tau, depth1, depth2, cross)
  } else {
    muhat = NF.no.DS(vehat, W[1 : n, ], W[-(1 : n), ], K, M, method2, depth2)
  }
  muhat.II = as.matrix(muhat[1 : n, ])
  
  ############################################################
  # our estimator and its variance 
  thetahatm1 = colMeans(muhat) + colMeans(vehat - muhat.II)
  thetahatm2 = colMeans(etahat) 
  thetahatm = thetahatm1 / thetahatm2
  thetahat = mean(thetahatm)
  
  psihatm1 = t( t( vehat - muhat.II ) / thetahatm2 )
  psihat1 = rowMeans(psihatm1)
  psihatm2 = t( t(muhat.II) / thetahatm2 - thetahatm1 * t(etahat.II) / thetahatm2 ^ 2 )
  psihat2 = rowMeans(psihatm2)
  sigmahatm = sqrt( colVar( psihatm1 + deltan * psihatm2 ) + deltan * (1 - deltan) * colVar(psihatm2) )
  sigmahat = sqrt( var(psihat1 + deltan * psihat2) + deltan * (1 - deltan) * var(psihat2) ) 
  ############################################################
  
  ############################################################
  # zero-imputation estimator and its variance 
  thetahatmzero1 = colMeans(vehat)
  thetahatmzero2 = colMeans(etahat)
  thetahatmzero = thetahatmzero1 / thetahatmzero2
  thetahatzero = mean(thetahatmzero)
  
  psihatmzero1 = t( t(vehat) / thetahatmzero2 )
  psihatzero1 = rowMeans(psihatmzero1)
  psihatmzero2 = - t( thetahatmzero1 * t(etahat.II) / thetahatmzero2 ^ 2 )
  psihatzero2 = rowMeans(psihatmzero2)
  sigmahatmzero = sqrt( colVar( psihatmzero1 + deltan * psihatmzero2 ) + deltan * (1 - deltan) * colVar(psihatmzero2) )
  sigmahatzero = sqrt( var( psihatzero1 + deltan * psihatzero2 ) + deltan * (1 - deltan) * var(psihatzero2) )
  ############################################################
  
  ############################################################
  # phase-II-only estimator and its variance
  thetahatOm1 = colMeans(vehat)
  thetahatOm2 = colMeans(etahat.II)
  thetahatOm = thetahatOm1 / thetahatOm2
  thetahatO = mean(thetahatOm)
  
  psihatOm = t( t(vehat) / thetahatOm2 - thetahatOm1 * t(etahat.II) / thetahatOm2 ^ 2 )
  psihatO = rowMeans(psihatOm)
  sigmahatOm = sqrt( colVar(psihatOm) )
  sigmahatO = sqrt( var(psihatO) )
  ############################################################
  
  result = list( "thetahatm" = thetahatm, "sdhatm" = sigmahatm / sqrt(n),
                 "thetahat" = thetahat,  "sdhat" = sigmahat / sqrt(n),
                 "thetahatmzero" = thetahatmzero, "sdhatmzero" = sigmahatmzero / sqrt(n),
                 "thetahatzero" = thetahatzero, "sdhatzero" = sigmahatzero / sqrt(n),
                 "thetahatOm" = thetahatOm, "sdhatOm" = sigmahatOm / sqrt(n),
                 "thetahatO" = thetahatO, "sdhatO" = sigmahatO / sqrt(n) )
  
  return(result)
  
}

# Conditional model fitting with cross-fitting
CM.fit = function(X, Y, K, M, method, tau = NULL, depth = NULL)
{
  
  n = length(Y)
  nK = floor(n / K)
  
  result = matrix(0, n, M)
  
  for(k in 1 : (K - 1))
  {
    
    ind = (k - 1) * nK + (1 : nK)
    
    Y.fit = Y[-ind]
    X.fit = as.matrix(X[-ind, ])
    X.predict = as.matrix(X[ind, ])
    
    result[ind, ] = CM.Est(Y.fit, X.fit, X.predict, method, tau, depth)
    
  }
  
  ind = ( ((K - 1) * nK) + 1 ) : n
  
  Y.fit = Y[-ind]
  X.fit = as.matrix(X[-ind, ])
  X.predict = as.matrix(X[ind, ])
  
  result[ind, ] = CM.Est(Y.fit, X.fit, X.predict, method, tau, depth)
  
  return(result)
  
}

# Nuisance function estimation without data-splitting
NF.no.DS = function(vehat, W.II, W.I, K, M, method2, depth2 = NULL)
{
  
  n = length(vehat)
  NO = nrow(W.I)
  nK = floor(n / K)
  
  result.I = matrix(0, NO, M)
  result.II = matrix(0, n, M)
  
  for(k in 1 : (K - 1))
  {
    
    ind = (k - 1) * nK + (1 : nK)
    which.fit = 1 : length(ind)
    
    vehat.fit = as.matrix(vehat[-ind, ])
    W.fit = W.II[-ind, ]
    W.predict = rbind(W.II[ind, ], W.I)
    
    NF.fit = apply( vehat.fit, 2, function(y) CM.Est(y, W.fit, W.predict, method2, tau = NULL, depth2) )
    
    result.I = result.I + NF.fit[-which.fit, ] / K
    result.II[ind, ] = NF.fit[which.fit, ]
    
  }
  
  ind = ( ((K - 1) * nK) + 1 ) : n
  which.fit = 1 : length(ind)
  
  vehat.fit = as.matrix(vehat[-ind, ])
  W.fit = W.II[-ind, ]
  W.predict = rbind(W.II[ind, ], W.I)
  
  NF.fit = apply( vehat.fit, 2, function(y) CM.Est(y, W.fit, W.predict, method2, tau = NULL, depth2) )
  
  result.I = result.I + NF.fit[-which.fit, ] / K
  result.II[ind, ] = NF.fit[which.fit, ]
  
  return( rbind(result.II, result.I) )
  
}

# Nuisance function estimation with data-splitting
NF.Est = function(Y.II, X.II, Z, W.I, K, M, method1, method2, loss, tau = NULL, depth1 = NULL, depth2 = NULL, cross)
{
  
  n = length(Y.II)
  NO = nrow(W.I)
  nK = floor(n / K)
  
  result.I = matrix(0, NO, M)
  result.II = matrix(0, n, M)
  
  W.II = cbind(as.matrix(Y.II), X.II)
  
  for(k in 1 : (K - 1))
  {
    
    ind = (k - 1) * nK + (1 : nK)
    which.fit = 1 : length(ind)
    
    Y.fit = Y.II[-ind]
    X.fit = as.matrix(X.II[-ind, ])
    Z.fit = as.matrix(Z[-ind, ])
    W.predict = rbind(W.II[ind, ], W.I)
    
    NF.fit = TS.NP.DS(Y.fit, X.fit, Z.fit, W.predict, M, method1, method2, loss, tau, depth1, depth2, cross)
    
    result.I = result.I + NF.fit[-which.fit, ] / K
    result.II[ind, ] = NF.fit[which.fit, ]
    
  }
  
  ind = ( ((K - 1) * nK) + 1 ) : n
  which.fit = 1 : length(ind)
  
  Y.fit = Y.II[-ind]
  X.fit = as.matrix(X.II[-ind, ])
  Z.fit = as.matrix(Z[-ind, ])
  W.predict = rbind(W.II[ind, ], W.I)
  
  NF.fit = TS.NP.DS(Y.fit, X.fit, Z.fit, W.predict, M, method1, method2, loss, tau, depth1, depth2, cross)
  
  result.I = result.I + NF.fit[-which.fit, ] / K
  result.II[ind, ] = NF.fit[which.fit, ]
  
  return( rbind(result.II, result.I) )
  
}

# Two-stage nonparametric estimator with data-splitting
# method1 is the method used for the first-stage nonparametric estimation
# method2 is the method used for the second-stage nonparametric estimation
TS.NP.DS = function(Y, X, Z, W.predict, M, method1, method2, loss, tau = NULL, depth1 = NULL, depth2 = NULL, cross)
{
  
  n = length(Y)
  
  n.half = floor(n / 2)
  
  W = cbind(as.matrix(Y), X)
  U = cbind(X, Z)
  
  result = matrix(0, nrow(W.predict), M)
  
  if(cross == T) 
  {
    J = 2
  } else {
    J = 1
  }
  
  for(j in 1 : J)
  {
    
    if(j == 1)
    {
      ind1 = 1 : n.half
      ind2 = (n.half + 1) : n
    } else {
      ind2 = 1 : n.half
      ind1 = (n.half + 1) : n
    }
    
    ghat.ds = CM.Est(Y[ind1], U[ind1, ], U[ind2, ], method1, tau, depth1)
    
    if(loss == "quadratic")
    {
      vehat.ds = quadratic(Y[ind2] - ghat.ds[, 1])
      result = result + CM.Est(vehat.ds, W[ind2, ], W.predict, method2, tau = NULL, depth2) / J
    }
    
    if(loss == "check")
    {
      
      vehat.ds = multiple.check(Y[ind2], ghat.ds, tau)
      result = result + apply( vehat.ds, 2, function(y) CM.Est(y, W[ind2, ], W.predict, method2, tau = NULL, depth2) ) / J
      
    }
    
  }
  
  return(result)
  
}

# Conditional model estimation
# depth is the interaction.depth paramter in the gbm function. 
# Setting it to 1 and 3 for {ghat, hhat} and muhat, respectively
CM.Est = function(Y.fit, X.fit, X.predict, method, tau = NULL, depth = NULL)
{
  
  if(method == "RF.mean")
  {
    
    model = randomForest(x = X.fit, y = Y.fit)
    result = predict(model, X.predict)
    
  }
  
  if(method == "LR.mean")
  {
    
    model = lm(Y.fit ~ X.fit)
    result = add.one(X.predict) %*% model $ coefficients
    
  }
  
  if(method == "GBM.mean")
  {
    
    data.fit = data.frame(Y.fit, X.fit)
    sink("aux")
    model = gbm(Y.fit ~ ., data = data.fit, var.monotone = numeric(ncol(X.fit)),
                distribution = "gaussian", n.trees = 1000, shrinkage = 0.01,
                interaction.depth = depth, bag.fraction = 0.5, train.fraction = 1,
                n.minobsinnode = 10, cv.folds = 5, keep.data = TRUE,
                verbose = F, n.cores = 1)
    best.iter = gbm.perf(model, plot.it = F, method = "cv")
    data.predict = data.frame( numeric(nrow(X.predict)), X.predict )
    result = predict(model, newdata = data.predict, n.trees = best.iter, type = "response")
    sink(NULL)
    
  }
  
  if(method == "GAM.mean")
  {
    
    p = ncol(X.fit)
    
    data.fit = data.frame(Y.fit, X.fit)
    names(data.fit) = c( "Y", paste0("X", 1 : p) )
    X.names = names(data.fit)[-1]
    gam.form = as.formula( paste0("Y ~ ", paste0("s(", X.names, ", bs = 'tp')", collapse = "+")) )
    model = mgcv::gam(gam.form, family = gaussian, data = data.fit,
                      method = "REML", optimizer = c("outer", "bfgs"))
    
    X.predict = data.frame(X.predict)
    names(X.predict) = X.names
    result = predict(model, X.predict, "response")
    
  }
  
  if(method == "XGB.mean")
  {
    
    sink("aux")
    model = xgboost(data = X.fit, label = Y.fit, max_depth = depth, eta = 0.3, nrounds = 1000, 
                    objective = "reg:squarederror", verbose = 0)
    sink(NULL)
    result = predict(model, X.predict)
  }
  
  if(method == "KS.mean")
  {
    
    sink("aux")
    h = npregbw(xdat = X.fit, ydat = Y.fit, nmulti = 1) $ bw
    result = npreg(txdat = X.fit, tydat = Y.fit, bws = h, exdat = X.predict) $ mean
    sink(NULL)
    
  }
  
  if(method == "RF.quantile")
  {
    
    model = quantregForest(x = X.fit, y = Y.fit)
    result = predict(model, X.predict, tau)
    
  }
  
  if(method == "LR.quantile")
  {
    
    model = rq(Y.fit ~ X.fit, tau = tau)
    result = add.one(X.predict) %*% model $ coefficients
    
  }
  
  if(method == "GAM.quantile")
  {
    
    p = ncol(X.fit)
    
    data.fit = data.frame(Y.fit, X.fit)
    names(data.fit) = c( "Y", paste0("X", 1 : p) )
    X.names = names(data.fit)[-1]
    gam.form = as.formula( paste0("Y ~ ", paste0("s(", X.names, ", bs = 'tp')", collapse = "+")) )
    sink("aux")
    model = qgam::mqgam( gam.form, data = data.fit, qu = tau,
                         argGam = list(method = "REML", optimizer = c("outer", "bfgs")) )
    sink(NULL)
    X.predict = data.frame(X.predict)
    names(X.predict) = X.names
    Y.predict = qdo(model, tau, predict, newdata = X.predict)
    result = matrix(unlist(Y.predict), ncol = length(tau))
    
  }
  
  return(as.matrix(result))
  
}
#########################################################################################################

outcome.generator = function(X, Z, model)
{
  
  N = nrow(X)
  
  xi = rnorm(N)
  
  if(model == "homo") 
  {
    
    Y = 2 * (Z[, 2]) ^ 3 / 5 + sqrt(2) * ( abs(Z[, 1]) + pnorm(X[, 1]) ) * xi 
    # Y = cos(X[, 1]) + 2 * (Z[, 2]) ^ 3 / 5 + sqrt(2) * ( abs(Z[, 1]) + pnorm(X[, 2]) ) * xi # good model
    
    
  }
  
  if(model == "hetero" )
  {
    
    Y = 2 * (X[, 1]) ^ 3 / 5 + sqrt(2) * ( abs(Z[, 1]) + pnorm(Z[, 2]) ) * xi
    # Y = cos(X[, 1]) + 2 * (X[, 2]) ^ 3 / 5 + sqrt(2) * ( abs(Z[, 1]) + pnorm(Z[, 2]) ) * xi # good model
    
  }
  
  return(Y)
  
}

# add a column of ones to the left of a matrix
add.one = function(X) { cbind(rep(1, nrow(X)), X) }

# calculate the length of an interval
CIL = function(x) x[2] - x[1]

# check whether theta0 is in the interval ( CI[1], CI[2] )
cover = function(theta0, CI) { (theta0 >= CI[1]) & (theta0 <= CI[2]) }

colVar = function(X) apply(X, 2, var)

# cqs stands for "fitted conditional quantiles"
multiple.check = function(Y, cqs, tau) { apply( rbind(tau, cqs), 2, function(x) check(Y - x[-1], x[1]) ) }

quadratic = function(x) { x ^ 2 }

check = function(x, tau) { x * ( tau - (x < 0) ) }

expit = function(x) { 1 / (1 + exp(-x)) }

bound = function(x) { c( max(0, x[1]), min(1, x[2]) ) }