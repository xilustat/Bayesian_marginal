rm(list=ls(all=TRUE))
library(mnormt)
library(SuppDists) # rinvGauss
library(miscTools) # colMeans
library(MASS)
library(MCMCpack)
source("LADBLSS.R")

# simulated data
n=200;m=40;
sig = matrix(0,m,m)

for (i in 1:m)
{
  for(j in 1:m)
  {
    sig[i,j] = 0.5^abs(i-j)
  }
}

set.seed(1234)
x = mvrnorm(n,rep(0,m),sig)
x = scale(x)
error = rnorm(n,0,1)
beta_true = rep(0,m)
beta_true[6:10] = runif(5,0.3,0.5)
beta_true[21:25] = runif(5,0.4,0.6)
beta_true[36:40] = runif(5,0.4,0.7)

y = x%*%beta_true + error 
y = as.matrix(y)

max.steps = 10000
betas1 <- LADBLSS(x, y, max.steps = 10000)
betas2 <- LADBLSS(x, y, max.steps = 10000)
betas3 <- LADBLSS(x, y, max.steps = 10000)
betas4 <- LADBLSS(x, y, max.steps = 10000)
betas5 <- LADBLSS(x, y, max.steps = 10000)


save(betas1,betas2,betas3,betas4,betas5,file = "PSRF_LADBLSS.RData")
