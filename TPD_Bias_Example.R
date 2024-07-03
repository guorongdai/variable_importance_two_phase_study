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

# method2 = "LR.mean"
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

r = 500

output = foreach(i = 1 : r, .combine = rbind) %dopar%
  {
    
    set.seed(789 * i)
    
    U = mvrnorm(N, numeric(p), Sigma)
    X = as.matrix( U[, 1 : p1] )
    Z = as.matrix( U[, -(1 : p1)] )
    Y = outcome.generator(X, Z, model)
    
    Y.center = ( Y > quantile(Y, 0.25) ) & ( Y < quantile(Y, 0.75) )
    pi = 0.3 * (1 - Y.center) + 0.1 * Y.center
    
    # pitilde = expit(X[, 1] / 2 + Y / 6 - 3 / 2)
    # pi = sapply(pitilde, function(x) min(x, 0.5))
    # pi = sapply(pi, function(x) max(x, 0.05))
    
    R = rbinom(N, 1, pi)
    
    Est = MVIM.TPD.B(Y, X, Z, R, pi, K, method1, method2, tau = tau, loss = loss, depth1, depth2)
    
    thetatilde = Est $ thetatilde 
    sdtilde = Est $ sdtilde
    CI = bound( thetatilde + c(-1, 1) * 1.96 * sdtilde )
    
    thetatildezero = Est $ thetatildezero
    sdtildezero = Est $ sdtildezero
    CIzero = bound( thetatildezero + c(-1, 1) * 1.96 * sdtildezero )
    
    thetatildeO = Est $ thetatildeO
    sdtildeO = Est $ sdtildeO
    CIO = bound( thetatildeO + c(-1, 1) * 1.96 * sdtildeO )
    
    all.CI = rbind(CI, CIzero, CIO)
    
    c( ( c(thetatilde, thetatildezero, thetatildeO) - theta0 ) ^ 2, 
       apply( all.CI, 1, CIL ), apply( all.CI, 1, function(x) cover(theta0, x) ) )
    
  }

registerDoSEQ()

output = apply(output, 2, mean)
output[1 : 3] = output[1 : 3] / output[1]

names(output) = c("MSE of thetatilde", "MSE of thetatildezero", "MSE of thetatildeO",
                  "CIL of thetatilde", "CIL of thetatildezero", "CIL of thetatildeO",
                  "CR of thetatilde", "CR of thetatildezero", "CR of thetatildeO")

output
