suppressMessages({
  
  if (!require(MASS)) install.packages('MASS', repos =
                                         'https://cran.revolutionanalytics.com/')
  if (!require(doParallel)) install.packages('parallel', repos =
                                               'https://cran.revolutionanalytics.com/')
  
  library(MASS)
  library(doParallel)
  
})

registerDoParallel(150) # set multicores
r = 1000
true.value = matrix(c(1, 0.9635653, 0.6385743, 0.7973178), 2, 2)
Y.quantile = matrix(c(-1.224684, -1.224849, 1.213393, 1.217346), 2, 2)

setwd("/home/guorongdai/TPD")
source("TPD_Functions.R")

N = 6000
p1 = 1
p2 = 2
p = p1 + p2

rho = 0.2
Sigma = matrix(0, p, p)
for(i in 1 : p) for(j in 1 : p) Sigma[i, j] = rho ^ abs(i - j)

model = "homo"

K = 10

method2 = "GBM.mean"

tau = (1 : 3) / 4

depth1 = 1 # interaction.depth parameter in the gbm function for estimating ghat and hhat
depth2 = 3 # interaction.depth parameter in the gbm function for estimating muhat
ds = T
cross = T

res.SRS = matrix(0, 24, 16) # results of the simple random sampling cases
res.BS = matrix(0, 24, 16) # results of the biased sampling cases

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
  
  for(j2 in 1 : 2) # 1 hetero; 2 homo
  {
    
    theta0 = true.value[j1, j2]
    
    if(j2 == 1) model = "hetero"
    if(j2 == 2) model = "homo"
    
    for(j3 in 1 : 2) # 1 SRS; 2 PS
    {
      
      for(j4 in 1 : 3) # phase-II design
      {
        
        if(j3 == 1)
        {
          
          n = c(300, 600, 1200)[j4]
          
          ss = 111
          
          output = foreach(i = 1 : r, .combine = rbind) %dopar%
            {
              
              set.seed(ss * i)
              
              U = mvrnorm(N, numeric(p), Sigma)
              X = as.matrix( U[, 1 : p1] )
              Z = as.matrix( U[, -(1 : p1)] )
              Y = outcome.generator(X, Z, model)
              
              Z.II = as.matrix( Z[1 : n, ] )
              
              MVIM.TPD(Y, X, Z.II, K, method1, method2, tau = tau, loss = loss, depth1, depth2, ds, cross)
              
            }
          
        }
        
        if(j3 == 2)
        {
          
          ss = 222
          
          output = foreach(i = 1 : r, .combine = rbind) %dopar%
            {
              
              set.seed(ss * i)
              
              U = mvrnorm(N, numeric(p), Sigma)
              X = as.matrix( U[, 1 : p1] )
              Z = as.matrix( U[, -(1 : p1)] )
              Y = outcome.generator(X, Z, model)
              
              if(j4 == 1)
              {
                # Y.center = ( Y > Y.quantile[j2, 1] ) & ( Y < Y.quantile[j2, 2] )
                Y.center = ( Y > quantile(Y, 0.25) ) & ( Y < quantile(Y, 0.75) )
                pi = 0.15 * (1 - Y.center) + 0.05 * Y.center
              }
              
              if(j4 == 2)
              {
                # Y.center = ( Y > Y.quantile[j2, 1] ) & ( Y < Y.quantile[j2, 2] )
                Y.center = ( Y > quantile(Y, 0.25) ) & ( Y < quantile(Y, 0.75) )
                pi = 0.3 * (1 - Y.center) + 0.1 * Y.center
              }
              
              if(j4 == 3)
              {
                pitilde = expit(X[, 1] / 2 + Y / 6 - 3 / 2)
                pi = sapply(pitilde, function(x) min(x, 0.5))
                pi = sapply(pi, function(x) max(x, 0.05))
              }
              
              R = rbinom(N, 1, pi)
              
              MVIM.TPD.B(Y, X, Z, R, pi, K, method1, method2, tau = tau, loss = loss, depth1, depth2)
              
            }
          
        }
        
        output.mean = apply(output, 2, mean)
        # output.mean[1 : 3] = output.mean[1 : 3] / output.mean[1] # calculate relative efficiencies
        output.sd = apply(output, 2, sd)
        output = rbind(output.mean, output.sd) * 100
        
        if(j3 == 1) res.SRS[ (j1 - 1) * 12 + (j2 - 1) * 6 + (j4 - 1) * 2 + (1 : 2), ] = output
        if(j3 == 2) res.BS[ (j1 - 1) * 12 + (j2 - 1) * 6 + (j4 - 1) * 2 + (1 : 2), ] = output
        
        print(c(j1, j2, j3, j4))
        
      }
      
    }
    
  }
  
}

registerDoSEQ()

print("###########SRS##########")
res.SRS
print("###########BS##########")
res.BS

write.csv(res.SRS, "res.SRS.csv", row.names = F)
write.csv(res.BS, "res.BS.csv", row.names = F)
