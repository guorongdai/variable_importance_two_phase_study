suppressMessages({
  
  if (!require(doParallel)) install.packages('parallel', repos =
                                               'https://cran.revolutionanalytics.com/')
  
  library(doParallel)
  
})

registerDoParallel(150) # set multicores
r = 1000
true.value = matrix(c(0.9992909, 0.9985872, 0.7131335, 0.8265324), 2, 2)

setwd("/home/guorongdai/TPD")
source("TPD_Functions.R")
dat = readRDS("dat.RDS")

Y = dat[, 1]

N = length(Y)

Y.quantile = quantile(Y, c(0.25, 0.75))

K = 10

method2 = "GBM.mean"

tau = (1 : 3) / 4

depth1 = 1 # interaction.depth parameter in the gbm function for estimating ghat and hhat
depth2 = 3 # interaction.depth parameter in the gbm function for estimating muhat
ds = T
cross = T

res.SRS = matrix(0, 12, 16) # results of the simple random sampling cases
res.BS = matrix(0, 12, 16) # results of the biased sampling cases

ss = 789

for(j1 in 1 : 2) # 1 mean; 2 quantiles
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
              
              MVIM.TPD(Y, X, Z.II, K, method1, method2, tau = tau, loss = loss, depth1, depth2, ds, cross)
              
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
              
              MVIM.TPD.B(Y, X, Z, R, pi, K, method1, method2, tau = tau, loss = loss, depth1, depth2)
              
            }
          
        }
        
        output.mean = apply(output, 2, mean)
        # output.mean[1 : 3] = output.mean[1 : 3] / output.mean[1] # calculate relative efficiencies
        output.sd = apply(output, 2, sd)
        output = rbind(output.mean, output.sd) * 100
        
        if(j3 == 1) res.SRS[ (j1 - 1) * 6 + (j4 - 1) * 2 + (1 : 2), ] = output
        if(j3 == 2) res.BS[ (j1 - 1) * 6 + (j4 - 1) * 2 + (1 : 2), ] = output
        
        print(c(j1, j3, j4))
        
      }
      
    }
    
  }
  
}

registerDoSEQ()

print("###########SRS##########")
res.SRS
print("###########BS##########")
res.BS

write.csv(res.SRS, "data.SRS.csv", row.names = F)
write.csv(res.BS, "data.BS.csv", row.names = F)
