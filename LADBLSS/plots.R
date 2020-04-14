
rm(list=ls(all=TRUE))

setwd("/Users/xilu/Desktop/LADBLSS")
load("PSRF_LADBLSS.RData")

library(mnormt)
library(VGAM) # rinv.gaussian
library(miscTools) # colMeans
library(MASS)
library(MCMCpack)
library(expm)
library(GIGrvg)
library(latex2exp)
#source("BLSS.R")


max.steps = 10000
d = 13  # the dth predictor
t1=as.matrix(betas1$t1[,d])
t1 = t1[seq(max.steps/2, max.steps,1),]
t2=as.matrix(betas2$t1[,d])
t2 = t2[seq(max.steps/2, max.steps,1),]
t3=as.matrix(betas3$t1[,d])
t3 = t3[seq(max.steps/2, max.steps,1),]
t4=as.matrix(betas4$t1[,d])
t4 = t4[seq(max.steps/2, max.steps,1),]
t5=as.matrix(betas5$t1[,d])
t5 = t5[seq(max.steps/2, max.steps,1),]
t11 = mcmc(t1)
t12 = mcmc(t2)
t13 = mcmc(t3)
t14 = mcmc(t4)
t15 = mcmc(t5)
xx = mcmc.list(t11,t12,t13,t14,t15)
gelman.diag(xx)
gelman.plot(xx,ylab = TeX('$\\beta_{1}$'), xlab = "BL-SS", ylim=c(1,1.5))

##################
# t1 = vSamples
##################

max.steps = 10000
t1=t2=t3=t4=t5=NULL
t11=t12=t13=t14=t15=NULL
d = 36  # the dth predictor
t1=as.matrix(betas1$t1[,d])
t1 = t1[seq(1, max.steps,1),]
t2=as.matrix(betas2$t1[,d])
t2 = t2[seq(1, max.steps,1),]
t3=as.matrix(betas3$t1[,d])
t3 = t3[seq(1, max.steps,1),]
t4=as.matrix(betas4$t1[,d])
t4 = t4[seq(1, max.steps,1),]
t5=as.matrix(betas5$t1[,d])
t5 = t5[seq(1, max.steps,1),]
t11 = mcmc(t1)
t12 = mcmc(t2)
t13 = mcmc(t3)
t14 = mcmc(t4)
t15 = mcmc(t5)
xx = mcmc.list(t11,t12,t13,t14,t15)
sum(abs(t1)+abs(t2)+abs(t3)+abs(t4)+abs(t5))
gelman.diag(xx)
gelman.plot(xx,ylab = "v",xlab = "LADBL", ylim=c(1,1.5))
temp = gelman.plot(xx,ylab = "v",xlab = "LADBL", ylim=c(1,1.5))
temp

##################
# t2 = sSamples
##################

max.steps = 10000
t1=t2=t3=t4=t5=NULL
t11=t12=t13=t14=t15=NULL
d = 36  # the dth predictor
t1=as.matrix(betas1$t2[,d])
t1 = t1[seq(1, max.steps,1),]
t2=as.matrix(betas2$t2[,d])
t2 = t2[seq(1, max.steps,1),]
t3=as.matrix(betas3$t2[,d])
t3 = t3[seq(1, max.steps,1),]
t4=as.matrix(betas4$t2[,d])
t4 = t4[seq(1, max.steps,1),]
t5=as.matrix(betas5$t2[,d])
t5 = t5[seq(1, max.steps,1),]
t11 = mcmc(t1)
t12 = mcmc(t2)
t13 = mcmc(t3)
t14 = mcmc(t4)
t15 = mcmc(t5)
xx = mcmc.list(t11,t12,t13,t14,t15)
sum(abs(t1)+abs(t2)+abs(t3)+abs(t4)+abs(t5))
gelman.diag(xx)
gelman.plot(xx,ylab = "s",xlab = "LADBL", ylim=c(1,1.5))
temp = gelman.plot(xx,ylab = "s", xlab = "LADBL", ylim=c(1,1.5))
temp

##################
# t3 = betaSamples
##################

max.steps = 10000
t1=t2=t3=t4=t5=NULL
t11=t12=t13=t14=t15=NULL
d = 41  # the dth predictor
t1=as.matrix(betas1$t3[,d])
t1 = t1[seq(1, max.steps,1),]
t2=as.matrix(betas2$t3[,d])
t2 = t2[seq(1, max.steps,1),]
t3=as.matrix(betas3$t3[,d])
t3 = t3[seq(1, max.steps,1),]
t4=as.matrix(betas4$t3[,d])
t4 = t4[seq(1, max.steps,1),]
t5=as.matrix(betas5$t3[,d])
t5 = t5[seq(1, max.steps,1),]
t11 = mcmc(t1)
t12 = mcmc(t2)
t13 = mcmc(t3)
t14 = mcmc(t4)
t15 = mcmc(t5)
xx = mcmc.list(t11,t12,t13,t14,t15)
sum(abs(t1)+abs(t2)+abs(t3)+abs(t4)+abs(t5))
gelman.diag(xx)
gelman.plot(xx,ylab = "beta",xlab = "LADBL", ylim=c(-100,100))
temp = gelman.plot(xx,ylab = TeX('$\\beta$'),xlab = "LADBL", ylim=c(1,1.5))
temp



##################
# t4 = tauSamples
##################

max.steps = 10000
t1=t2=t3=t4=t5=NULL
t11=t12=t13=t14=t15=NULL
d = 1  # the dth predictor
t1=as.matrix(betas1$t4[,d])
t1 = t1[seq(1, max.steps,1),]
t2=as.matrix(betas2$t4[,d])
t2 = t2[seq(1, max.steps,1),]
t3=as.matrix(betas3$t4[,d])
t3 = t3[seq(1, max.steps,1),]
t4=as.matrix(betas4$t4[,d])
t4 = t4[seq(1, max.steps,1),]
t5=as.matrix(betas5$t4[,d])
t5 = t5[seq(1, max.steps,1),]
t11 = mcmc(t1)
t12 = mcmc(t2)
t13 = mcmc(t3)
t14 = mcmc(t4)
t15 = mcmc(t5)
xx = mcmc.list(t11,t12,t13,t14,t15)
sum(abs(t1)+abs(t2)+abs(t3)+abs(t4)+abs(t5))
gelman.diag(xx)
temp=gelman.plot(xx,ylab = TeX('$\\tau$'),xlab = "LADBL", ylim=c(1,1.5))
temp


# gelman.plot(xx,ylab = TeX('$\\beta_{1}$'),xlab = "BL-SS")


##################
# t5 = eta2Samples
##################

max.steps = 10000
t1=t2=t3=t4=t5=NULL
t11=t12=t13=t14=t15=NULL
temp=xx=NULL
d = 1 # the dth predictor
t1=as.matrix(betas1$t5[,d])
t1 = t1[seq(1, max.steps,1),]
t2=as.matrix(betas2$t5[,d])
t2 = t2[seq(1, max.steps,1),]
t3=as.matrix(betas3$t5[,d])
t3 = t3[seq(1, max.steps,1),]
t4=as.matrix(betas4$t5[,d])
t4 = t4[seq(1, max.steps,1),]
t5=as.matrix(betas5$t5[,d])
t5 = t5[seq(1, max.steps,1),]
t11 = mcmc(t1)
t12 = mcmc(t2)
t13 = mcmc(t3)
t14 = mcmc(t4)
t15 = mcmc(t5)
xx = mcmc.list(t11,t12,t13,t14,t15)
sum(abs(t1)+abs(t2)+abs(t3)+abs(t4)+abs(t5))
gelman.diag(xx)
temp=gelman.plot(xx,ylab = TeX('$\\eta^{2}$'),xlab = "LADBL")
temp


##################
# The Plot
# t3 = betaSamples
##################

for(i in 1:100){
  
  t1=t2=t3=t4=t5=NULL
  t11=t12=t13=t14=t15=NULL
  d = i  # the dth predictor
  t1=as.matrix(betas1$t3[,d])
  t1 = t1[seq(1, max.steps,1),]
  t2=as.matrix(betas2$t3[,d])
  t2 = t2[seq(1, max.steps,1),]
  t3=as.matrix(betas3$t3[,d])
  t3 = t3[seq(1, max.steps,1),]
  t4=as.matrix(betas4$t3[,d])
  t4 = t4[seq(1, max.steps,1),]
  t5=as.matrix(betas5$t3[,d])
  t5 = t5[seq(1, max.steps,1),]
  t11 = mcmc(t1)
  t12 = mcmc(t2)
  t13 = mcmc(t3)
  t14 = mcmc(t4)
  t15 = mcmc(t5)
  xx = mcmc.list(t11,t12,t13,t14,t15)
  sum(abs(t1)+abs(t2)+abs(t3)+abs(t4)+abs(t5))
  gelman.diag(xx)
  gelman.plot(xx,main=paste("The",d,"th"),ylab = TeX('$\\beta$'),xlab = "LADBL", ylim=c(1,1.5))
  temp = gelman.plot(xx,main=paste("The",d,"th"),ylab = TeX('$\\beta$'),xlab = "LADBL", ylim=c(1,1.5))
  temp
  
  Sys.sleep(0.3)
}


##################
# The Plot
# t1 = vSamples
##################

for(i in 1:100){
  
  t1=t2=t3=t4=t5=NULL
  t11=t12=t13=t14=t15=NULL
  d = i  # the dth predictor
  t1=as.matrix(betas1$t1[,d])
  t1 = t1[seq(1, max.steps,1),]
  t2=as.matrix(betas2$t1[,d])
  t2 = t2[seq(1, max.steps,1),]
  t3=as.matrix(betas3$t1[,d])
  t3 = t3[seq(1, max.steps,1),]
  t4=as.matrix(betas4$t1[,d])
  t4 = t4[seq(1, max.steps,1),]
  t5=as.matrix(betas5$t1[,d])
  t5 = t5[seq(1, max.steps,1),]
  t11 = mcmc(t1)
  t12 = mcmc(t2)
  t13 = mcmc(t3)
  t14 = mcmc(t4)
  t15 = mcmc(t5)
  xx = mcmc.list(t11,t12,t13,t14,t15)
  sum(abs(t1)+abs(t2)+abs(t3)+abs(t4)+abs(t5))
  gelman.diag(xx)
  gelman.plot(xx,main=paste("The",d,"th"),ylab = "v",xlab = "LADBL", ylim=c(1,1.5))
  temp = gelman.plot(xx,main=paste("The",d,"th"),ylab = "v",xlab = "LADBL", ylim=c(1,1.5))
  temp
  
  Sys.sleep(0.3)
}

##################
# The Plot
# t2 = sSamples
##################


for(i in 1:100){
  
  t1=t2=t3=t4=t5=NULL
  t11=t12=t13=t14=t15=NULL
  d = i  # the dth predictor
  t1=as.matrix(betas1$t2[,d])
  t1 = t1[seq(1, max.steps,1),]
  t2=as.matrix(betas2$t2[,d])
  t2 = t2[seq(1, max.steps,1),]
  t3=as.matrix(betas3$t2[,d])
  t3 = t3[seq(1, max.steps,1),]
  t4=as.matrix(betas4$t2[,d])
  t4 = t4[seq(1, max.steps,1),]
  t5=as.matrix(betas5$t2[,d])
  t5 = t5[seq(1, max.steps,1),]
  t11 = mcmc(t1)
  t12 = mcmc(t2)
  t13 = mcmc(t3)
  t14 = mcmc(t4)
  t15 = mcmc(t5)
  xx = mcmc.list(t11,t12,t13,t14,t15)
  sum(abs(t1)+abs(t2)+abs(t3)+abs(t4)+abs(t5))
  gelman.diag(xx)
  gelman.plot(xx,main=paste("The",d,"th"),ylab = "s",xlab = "LADBL", ylim=c(1,1.5))
  temp = gelman.plot(xx,main=paste("The",d,"th"),ylab = "s",xlab = "LADBL", ylim=c(1,1.5))
  temp
  
  Sys.sleep(0.3)
}


