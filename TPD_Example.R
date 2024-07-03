suppressMessages({
  
  if (!require(MASS)) install.packages('MASS', repos =
                                         'https://cran.revolutionanalytics.com/')
  if (!require(doParallel)) install.packages('parallel', repos =
                                               'https://cran.revolutionanalytics.com/')
  
  library(MASS)
  library(doParallel)
  
})

registerDoParallel(12) # set multicores

source("/Users/daiguorong/Library/CloudStorage/Dropbox/TPD/Simulation/TPD_Functions.R")

N = 6000
n = 1200 # n = 300, 600 or 1200
p1 = 1
p2 = 2
p = p1 + p2

rho = 0.2
Sigma = matrix(0, p, p)
for(i in 1 : p) for(j in 1 : p) Sigma[i, j] = rho ^ abs(i - j)

model = "homo"

K = 10

method1 = "GAM.mean"
# method1 = "GAM.quantile"

if(method1 == "GAM.mean") loss = "quadratic"
if(method1 == "GAM.quantile") loss = "check"

method2 = "GBM.mean"

tau = (1 : 3) / 4

depth1 = 1 # interaction.depth parameter in the gbm function for estimating ghat and hhat
depth2 = 3 # interaction.depth parameter in the gbm function for estimating muhat
ds = T
cross = T

######################
# calculate theta0 for the quantile cases
# tau = c(0.25, 0.50, 0.75)
# theta0 = 0.9635653 # thetam = c(0.9476285, 1, 0.9430347) # hetero
# theta0 = 0.7973178 # thetam = c(0.7856039, 0.8238681, 0.7824814) # homo
######################

######################
# calculate theta0 for the mean cases
# theta0 = 1
theta0 = 0.6385743 # p1 = 1, homo
######################

r = 1000

output = foreach(i = 1 : r, .combine = rbind) %dopar%
  {
    
    set.seed(777 * i) # others
    # set.seed(789 * i) # mean: hetero 600 & 1200; quantiles: hetero 1200
    
    U = mvrnorm(N, numeric(p), Sigma)
    X = as.matrix( U[, 1 : p1] )
    Z = as.matrix( U[, -(1 : p1)] )
    Y = outcome.generator(X, Z, model)
    
    Z.II = as.matrix( Z[1 : n, ] )
    
    Est = MVIM.TPD(Y, X, Z.II, K, method1, method2, tau = tau, loss = loss, depth1, depth2, ds, cross)
    
    thetahat = Est $ thetahat 
    sdhat = Est $ sdhat
    CI = bound(thetahat + c(-1, 1) * 1.96 * sdhat)
    
    thetahatzero = Est $ thetahatzero
    sdhatzero = Est $ sdhatzero
    CIzero = bound(thetahatzero + c(-1, 1) * 1.96 * sdhatzero)
    
    thetahatO = Est $ thetahatO
    sdhatO = Est $ sdhatO
    CIO = bound(thetahatO + c(-1, 1) * 1.96 * sdhatO)
    
    c( ( c(thetahat, thetahatzero, thetahatO) - theta0 ) ^ 2,
       apply( rbind(CI, CIzero, CIO), 1, CIL ), apply( rbind(CI, CIzero, CIO), 1, function(x) cover(theta0, x) ) )
    
  }

registerDoSEQ()

output = apply(output, 2, mean)
output[1 : 3] = output[1 : 3] / output[1]

names(output) = c("MSE of thetahat", "MSE of thetahatzero", "MSE of thetahatO",
                  "CIL of thetahat", "CIL of thetahatzero", "CIL of thetahatO",
                  "CR of thetahat", "CR of thetahatzero", "CR of thetahatO")

output