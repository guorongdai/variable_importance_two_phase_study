suppressMessages({
  
  if (!require(doParallel)) install.packages('parallel', repos =
                                               'https://cran.revolutionanalytics.com/')
  
  library(doParallel)
  
})

registerDoParallel(255) # set multicores
r = 500
true.value = matrix(c(0.9992909, 0.9985872, 0.7131335, 0.8265324), 2, 2)

source("/home/guorongdai/TPD/TPD_Functions.R")
dat = readRDS("/home/guorongdai/TPD/Data_Analysis/dat.RDS")

Y = dat[, 1]

N = length(Y)

Y.quantile = quantile(Y, c(0.25, 0.75))

K = 10

method2 = "GBM.mean"

tau = (1 : 3) / 4

depth2 = 3 # interaction.depth parameter in the gbm function for estimating muhat
ds = T
cross = T

RE = matrix(0, 6, 8)
CI = matrix(0, 12, 12)

ss = 789

for(j1 in 2) # 1 mean; 2 quantiles
{
  
  if(j1 == 1)
  {
    method1 = "GAM.mean"
    loss = "quadratic"
  } else {
    method1 = "GAM.quantile"
    loss = "check"
  }
  
  for(j2 in 1) # 1: scenario 1, X is demographic and Z is dietary; 2: scenario 2, the other way around
  {
    
    theta0 = true.value[j1, j2]
    
    X = dat[, 2 : 3]
    Z = as.matrix(dat[, 4 : 5])
    
    if(j2 == 2)
    {
      A = X
      X = Z
      Z = A
    }
    
    for(j3 in 1 : 2) # 1 SRS; 2 PS
    {
      
      for(j4 in 1 : 3) # phase-II design
      {
        
        if(j3 == 1)
        {
          
          n = c(300, 600, 1200)[j4]
          
          output = foreach(i = 1 : r, .combine = rbind) %dopar%
            {
              
              set.seed(ss * i)
              
              II.no = sample(1 : N, n)
              Y = c(Y[II.no], Y[-II.no])
              X = as.matrix( rbind(X[II.no, ], X[-II.no, ]) )
              Z.II = as.matrix( Z[II.no, ] )
              
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
                 apply( rbind(CI, CIzero, CIO), 1, CIL ), apply( rbind(CIzero, CI, CIO), 1, function(x) cover(theta0, x) ) )
              
            }
          
        }
        
        if(j3 == 2)
        {
          
          output = foreach(i = 1 : r, .combine = rbind) %dopar%
            {
              
              if(j4 == 1)
              {
                Y.center = ( Y > Y.quantile[1] ) & ( Y < Y.quantile[2] )
                pi = 0.15 * (1 - Y.center) + 0.05 * Y.center
              }
              
              if(j4 == 2)
              {
                Y.center = ( Y > Y.quantile[1] ) & ( Y < Y.quantile[2] )
                pi = 0.3 * (1 - Y.center) + 0.1 * Y.center
              }
              
              if(j4 == 3)
              {
                pitilde = expit( (X[, 1] + X[, 2]) / 3 + Y / 6 - 3 / 2 )
                pi = sapply(pitilde, function(x) min(x, 0.5))
                pi = sapply(pi, function(x) max(x, 0.05))
              }
              
              set.seed(ss * i)
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
              
              c( ( c(thetatilde, thetatildezero, thetatildeO) - theta0 ) ^ 2,
                 apply( rbind(CI, CIzero, CIO), 1, CIL ), apply( rbind(CIzero, CI, CIO), 1, function(x) cover(theta0, x) ) )
              
            }
          
        }
        
        output = apply(output, 2, mean)
        output[1 : 3] = output[1 : 3] / output[1]
        
        RE[ (j3 - 1) * 3 + j4, (j1 - 1) * 4 + (j2 - 1) * 2 + (1 : 2) ] = output[2 : 3]
        CI[ (j1 - 1) * 6 + (j3 - 1) * 3 + j4, (j2 - 1) * 6 + (1 : 6) ] = output[4 : 9] * 100
        
        print(c(j1, j3, j4))
        
      }
      
    }
    
  }
  
}

registerDoSEQ()

print("###########RE##########")
RE
print("###########CI##########")
CI

write.csv(RE, "RE.csv", row.names = F)
write.csv(CI, "CI.csv", row.names = F)
