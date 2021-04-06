library(expm)
library(latex2exp)

auc<-function(x,y)
{
  n=length(x)
  a=rep(0,(n-1))
  for (i in 1:(n-1)) 
  {
    a[i]=(abs(x[i+1]-x[i]))*((y[i+1]+y[i])/2)
  }
  return(sum(a))
}

auc<-function(x,y)
{
  n=length(x)
  b=0
  for (i in 1:(n-1)) 
  {
    a=(abs(x[i+1]-x[i]))*((y[i+1]+y[i])/2)
    b=b+a
  }
  return(b)
}

#BL
load("roc_BL.RData")
fpr1=fpr_BL; tpr1=tpr_BL; auc1=auc_BL
fpr2=fpr_BL; tpr2=tpr_BL; auc2=auc_BL
fpr3=fpr_BL; tpr3=tpr_BL; auc3=auc_BL
fpr_BL=rbind(fpr1,fpr2,fpr3); tpr_BL=rbind(tpr1,tpr2,tpr3); auc_BL=c(auc1,auc2,auc3)

#auc_BL=c(auc1,auc2,auc3,auc4,auc5,auc6)
#auc_BL=c(auc1,auc2,auc3,auc4,auc5,auc6,auc7,auc8,auc9,auc10)

fpr=colMeans(fpr_BL); tpr=colMeans(tpr_BL); 
fpr_BL=fpr; tpr_BL=tpr;
#auc_BL 
mean(auc_BL)
sd(auc_BL) 

#fpr_BL=c(fpr_BL,fpr); 
#tpr_BL=c(tpr_BL,tpr);

matrix=cbind(fpr_BL,tpr_BL)
matrix = matrix[order(matrix[,1]),]
m1=matrix
fpr_BL=m1[,1]; tpr_BL=m1[,2]


save(tpr_BL, fpr_BL, auc_BL, file = "roc_BL.RData")
auc(fpr_BL,tpr_BL)
plot(fpr_BL,tpr_BL, type = "l",lty=2)

#BLSS
load("roc_BLSS.RData")
fpr1=fpr_BLSS; tpr1=tpr_BLSS; auc1=auc_BLSS
fpr2=fpr_BLSS; tpr2=tpr_BLSS; auc2=auc_BLSS

fpr_BLSS=rbind(fpr1,fpr2); tpr_BLSS=rbind(tpr1,tpr2); auc_BLSS=c(auc1,auc2)

fpr=colMeans(fpr_BLSS); tpr=colMeans(tpr_BLSS); 
fpr_BLSS=fpr; tpr_BLSS=tpr;
#auc_BLSS 
mean(auc_BLSS)
sd(auc_BLSS) 


matrix2=cbind(fpr_BLSS,tpr_BLSS)
matrix2 = matrix2[order(matrix2[,1]),]
m2=matrix2
fpr_BLSS=m2[,1]; tpr_BLSS=m2[,2]

save(tpr_BLSS, fpr_BLSS,auc_BLSS, file = "roc_BLSS.RData")
plot(fpr_BLSS,tpr_BLSS,type = "l",lty=2)
auc(fpr_BLSS,tpr_BLSS)

#LADBL
load("roc_LADBL.RData")
fpr1=fpr_LADBL; tpr1=tpr_LADBL; auc1=auc_LADBL
fpr2=fpr_LADBL; tpr2=tpr_LADBL; auc2=auc_LADBL


fpr_LADBL=rbind(fpr1,fpr2)
tpr_LADBL=rbind(tpr1,tpr2)
auc_LADBL=c(auc1,auc2)

fpr=colMeans(fpr_LADBL); tpr=colMeans(tpr_LADBL); 
fpr_LADBL=fpr; tpr_LADBL=tpr;
#auc_LADBL
mean(auc_LADBL)
sd(auc_LADBL) 


matrix3=cbind(fpr_LADBL,tpr_LADBL)
matrix3 = matrix3[order(matrix3[,1]),]
m3=matrix3
fpr_LADBL=m3[,1]; tpr_LADBL=m3[,2]

save(tpr_LADBL, fpr_LADBL, auc_LADBL, file = "roc_LADBL.RData")
auc(fpr_LADBL,tpr_LADBL)
plot(fpr_LADBL,tpr_LADBL,type = "l",lty=2)


#LADBLSS
load("roc_LADBLSS.RData")
fpr1=fpr_LADBLSS; tpr1=tpr_LADBLSS; auc1=auc_LADBLSS
fpr2=fpr_LADBLSS; tpr2=tpr_LADBLSS; auc2=auc_LADBLSS
#fpr3=fpr_LADBLSS; tpr3=tpr_LADBLSS; auc3=auc_LADBLSS


fpr_LADBLSS=rbind(fpr1,fpr2)
tpr_LADBLSS=rbind(tpr1,tpr2)
auc_LADBLSS=c(auc1,auc2)


fpr=colMeans(fpr_LADBLSS); tpr=colMeans(tpr_LADBLSS); 
fpr_LADBLSS=fpr; tpr_LADBLSS=tpr;
#auc_LADBLSS
mean(auc_LADBLSS)
sd(auc_LADBLSS)


matrix4=cbind(fpr_LADBLSS,tpr_LADBLSS)
matrix4 = matrix4[order(matrix4[,1]),]
m4=matrix4[,];
fpr_LADBLSS=m4[,1]; tpr_LADBLSS=m4[,2]

save(tpr_LADBLSS, fpr_LADBLSS, auc_LADBLSS, file = "roc_LADBLSS.RData")
auc(fpr_LADBLSS,tpr_LADBLSS)
plot(fpr_LADBLSS,tpr_LADBLSS,type = "l",lty=2)


library(expm)
library(latex2exp)
setwd("/Users/xilu/Desktop/simulation 2/ROC")
load("roc_BL.RData"); load("roc_BLSS.RData"); load("roc_LADBL.RData"); load("roc_LADBLSS.RData"); 
par(mfrow=c(1,1))
plot(fpr_BL,tpr_BL,type = "l", col="red",lty=2,xlim = c(0,1),ylim = c(0,1), xlab = "FPR",ylab = "TPR", main=TeX('error=80%n(0,1)+20%Cauchy(0,1)'))
lines(fpr_BLSS,tpr_BLSS,type = "l", col="blue",lty=2)
lines(fpr_LADBL,tpr_LADBL,type = "l", col="green",lty=2)
lines(fpr_LADBLSS,tpr_LADBLSS,type = "l", col="black",lty=2)

legend("bottomright", lty =c(2,2,2,2), col = c("red","blue","green","black"), legend = c("BL","BLSS","LADBL","LADBLSS"),cex = 0.5)
auc(fpr_BL,tpr_BL) #0.8593902
auc(fpr_BLSS,tpr_BLSS) #0.8771207
auc(fpr_LADBL,tpr_LADBL) #0.8581815
auc(fpr_LADBLSS,tpr_LADBLSS) #0.8639152

legend("bottom",legend = c("AUC","BL=0.7636","BLSS=0.7777","LADBL=0.7991","LADBLSS=0.8166"),cex = 0.5)

TeX('error=10%laplace(0,1)+90%laplace(0,$\\sqrt{5}$)')














