library(roben)


iter=10000
TP_G =c(); FP_G =c(); TP_GXE =c(); FP_GXE =c(); 

for (i in 1:30) 
{
  
  coef_G = c(); coef_GXE = c()
  
  for(j in 1:p)
  {
    X=g[,j]; Y=y; E=e; clin = c;
    fit=roben(X, Y, E, clin, iterations = iter)
    coef_G=c(coef_G,fit$coefficient$GE[1]) 
    coef_GXE = c(coef_GXE,fit$coefficient$GE[2:5])
  }
  
  TP_G = c(TP_G,tp(beta_true,coef_G)); FP_G = c(FP_G,fp(beta_true,coef_G)); 
  TP_GXE = c(TP_GXE,tp(eta_true,coef_GXE)); FP_GXE = c(FP_GXE,fp(eta_true,coef_GXE));   
}



sel=GxESelection(fit)
pos=which(sel$indicator != 0)
tp=length(intersect(which(coeff$GE != 0), pos))
fp=length(pos)-tp
list(tp=tp, fp=fp)


