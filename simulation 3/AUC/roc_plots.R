library(expm)
library(latex2exp)

setwd("/Users/xilu/Desktop/marginal_bayes/simulation/simulation 3/AUC/data2")
load("roc_BL.RData"); load("roc_BLSS.RData"); load("roc_LADBL.RData"); load("roc_LADBLSS.RData"); 

fpr_BL_1 = fpr_BL; tpr_BL_1 = tpr_BL; fpr_BLSS_1 = fpr_BLSS; tpr_BLSS_1 = tpr_BLSS; 
fpr_LADBL_1 = fpr_LADBL; tpr_LADBL_1 = tpr_LADBL; fpr_LADBLSS_1 = fpr_LADBLSS; tpr_LADBLSS_1 = tpr_LADBLSS; 

setwd("/Users/xilu/Desktop/marginal_bayes/simulation/simulation 3/AUC/data3")
load("roc_BL.RData"); load("roc_BLSS.RData"); load("roc_LADBL.RData"); load("roc_LADBLSS.RData"); 

fpr_BL_2 = fpr_BL; tpr_BL_2 = tpr_BL; fpr_BLSS_2 = fpr_BLSS; tpr_BLSS_2 = tpr_BLSS; 
fpr_LADBL_2 = fpr_LADBL; tpr_LADBL_2 = tpr_LADBL; fpr_LADBLSS_2 = fpr_LADBLSS; tpr_LADBLSS_2 = tpr_LADBLSS; 


setwd("/Users/xilu/Desktop/marginal_bayes/simulation/simulation 3/AUC/data4")
load("roc_BL.RData"); load("roc_BLSS.RData"); load("roc_LADBL.RData"); load("roc_LADBLSS.RData"); 

fpr_BL_3 = fpr_BL; tpr_BL_3 = tpr_BL; fpr_BLSS_3 = fpr_BLSS; tpr_BLSS_3 = tpr_BLSS; 
fpr_LADBL_3 = fpr_LADBL; tpr_LADBL_3 = tpr_LADBL; fpr_LADBLSS_3 = fpr_LADBLSS; tpr_LADBLSS_3 = tpr_LADBLSS; 

setwd("/Users/xilu/Desktop/marginal_bayes/simulation/simulation 3/AUC/data5")
load("roc_BL.RData"); load("roc_BLSS.RData"); load("roc_LADBL.RData"); load("roc_LADBLSS.RData"); 

fpr_BL_4 = fpr_BL; tpr_BL_4 = tpr_BL; fpr_BLSS_4 = fpr_BLSS; tpr_BLSS_4 = tpr_BLSS; 
fpr_LADBL_4 = fpr_LADBL; tpr_LADBL_4 = tpr_LADBL; fpr_LADBLSS_4 = fpr_LADBLSS; tpr_LADBLSS_4 = tpr_LADBLSS; 

par(mfrow=c(2,2))
plot(fpr_BL_1,tpr_BL_1,type = "l", col="red",lty=2,xlim = c(0,1),ylim = c(0,1), xlab = "FPR",ylab = "TPR", main=TeX('error=t(2)'))
lines(fpr_BLSS_1,tpr_BLSS_1,type = "l", col="blue",lty=2)
lines(fpr_LADBL_1,tpr_LADBL_1,type = "l", col="green",lty=2)
lines(fpr_LADBLSS_1,tpr_LADBLSS_1,type = "l", col="black",lty=2)

legend("bottomright", lty =c(2,2,2,2), col = c("red","blue","green","black"), legend = c("BL","BLSS","LADBL","LADBLSS"),cex = 0.5)


plot(fpr_BL_2,tpr_BL_2,type = "l", col="red",lty=2,xlim = c(0,1),ylim = c(0,1), xlab = "FPR",ylab = "TPR", main=TeX('error=Lognormal(0,2)'))
lines(fpr_BLSS_2,tpr_BLSS_2,type = "l", col="blue",lty=2)
lines(fpr_LADBL_2,tpr_LADBL_2,type = "l", col="green",lty=2)
lines(fpr_LADBLSS_2,tpr_LADBLSS_2,type = "l", col="black",lty=2)

legend("bottomright", lty =c(2,2,2,2), col = c("red","blue","green","black"), legend = c("BL","BLSS","LADBL","LADBLSS"),cex = 0.5)

plot(fpr_BL_3,tpr_BL_3,type = "l", col="red",lty=2,xlim = c(0,1),ylim = c(0,1), xlab = "FPR",ylab = "TPR", main=TeX('error=90%n(0,1)+10%Cauchy(0,1)'))
lines(fpr_BLSS_3,tpr_BLSS_3,type = "l", col="blue",lty=2)
lines(fpr_LADBL_3,tpr_LADBL_3,type = "l", col="green",lty=2)
lines(fpr_LADBLSS_3,tpr_LADBLSS_3,type = "l", col="black",lty=2)

legend("bottomright", lty =c(2,2,2,2), col = c("red","blue","green","black"), legend = c("BL","BLSS","LADBL","LADBLSS"),cex = 0.5)

plot(fpr_BL_4,tpr_BL_4,type = "l", col="red",lty=2,xlim = c(0,1),ylim = c(0,1), xlab = "FPR",ylab = "TPR", main=TeX('error=80%n(0,1)+20%Cauchy(0,1)'))
lines(fpr_BLSS_4,tpr_BLSS_4,type = "l", col="blue",lty=2)
lines(fpr_LADBL_4,tpr_LADBL_4,type = "l", col="green",lty=2)
lines(fpr_LADBLSS_4,tpr_LADBLSS_4,type = "l", col="black",lty=2)

legend("bottomright", lty =c(2,2,2,2), col = c("red","blue","green","black"), legend = c("BL","BLSS","LADBL","LADBLSS"),cex = 0.5)






