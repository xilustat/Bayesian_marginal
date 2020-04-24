LADBLSS = function(x, y, max.steps)
{
  n = nrow(x)
  p = ncol(x)
  vsample = matrix(0,max.steps,n)
  ssample = matrix(0,max.steps,p)
  betasample = matrix(0,max.steps,p)
  tausample = matrix(0,max.steps,1)
  eta2sample = matrix(0,max.steps,1)
  pisample = matrix(0,max.steps,1)
  SS = matrix(0,max.steps,p)
  v = rep(1,n)
  s = rep(1,p)
  beta = rep(1,p)
  z = rep(1,p)
  sg = rep(1,p)
  tau = 1
  eta2 = 1
  theta=0.5
  e2 = 2/(theta*(1-theta))
  a = 1e-1
  b = 1e-1
  c = 1e-1
  d = 1e-1
  pi=1/2
  
  for(k in 1:max.steps)
{
    iter = 0
    
    #sample v
    #res = y-x%*%beta
    #lambda = 2*tau
    #mu = sqrt(lambda*e2/(tau*res^2))
    #index = 1:n
    #flag=1
    #while(flag){
      #inv_v= rinvGauss(length(index),nu=mu[index],lambda=lambda)
      #flag = any(inv_v<=0)|any(is.na(inv_v))|any(is.infinite(inv_v))
      #v[index[inv_v>0]] = 1/inv_v[inv_v>0]
      #index = setdiff(index,index[inv_v>0])
    #}
    #vsample[k,] = v
    
    #sample v
    for(i in 1:n)
    {
    res = y-x%*%beta
    lambda = 2*tau
    mu = sqrt(lambda*e2/(tau*res[i]^2))
    flag=1
    while(flag){
      inv_v= rinvGauss(1,nu=mu,lambda=lambda)
      flag = any(inv_v<=0)|any(is.na(inv_v))|any(is.infinite(inv_v))
    }
    v[i]=1/inv_v
    }
    vsample[k,] = v
    
    iter=iter+1
    
    #sample beta
    for(j in 1:p){
      A = x[,j]^2/v
      invsigma2 = tau*sum(A)/e2+1/s[j]
      sigma2 = 1/invsigma2
      y_j = as.vector(y -x[,-j]%*%beta[-j])
      B = y_j*x[,j]/v
      BB = tau*sum(B)/e2
      mu = BB*sigma2
      c = (1/2)*(sqrt(sigma2)*BB)^2
      d = s[j]^(-1/2)*sqrt(sigma2)*exp(c)
      l = pi/(pi+(1-pi)*d)
      u = runif(1)
      if(u<l){
        beta[j] = 0; z[j]=0; sg[j]=1
      }
      else{
        beta[j] = rnorm(1,mean=mu,sd=sqrt(sigma2));z[j]=1;sg[j]=0
      }
    }
    betasample[k,]=beta
     
    iter=iter+1
    
    #sample s
    for(j in 1:p){
      if(beta[j]==0){
        s[j] = rexp(1,rate=eta2/2)
      }
      else{
        lambda_1 = eta2
        mu_1 = sqrt(lambda_1/beta[j]^2)
        flag=1
        while(flag){
          invs = rinvGauss(1,nu=mu_1,lambda=lambda_1)
          flag = (invs<=0)|(is.na(invs))|(is.infinite(invs))
        }
        s[j]=1/invs
      }
    }
    ssample[k,] = s
    
    iter=iter+1
    
    #sample tau
    res = y-x%*%beta
    vec = res^2/(2*e2*v)+v
    rate = sum(vec)+b
    shape = a+3*n/2
    tau = rgamma(1,shape=shape,rate=rate)
    tausample[k,] = tau
    
    iter=iter+1
    
    #sample eta2
    rate_1 = sum(s/2)+d
    shape_1 = p+c
    eta2 = rgamma(1,shape=shape_1,rate=rate_1)
    eta2sample[k,] = eta2
    
    iter=iter+1
    
    #sample pi
    shape1 = c + p - sum(z)
    shape2 = d + sum(z)
    pi = rbeta(1,shape1=shape1, shape2=shape2)
    pisample[k,] = pi
    
    iter=iter+1
    
    #spike and slab
    SS[k,] = sg
  }
  t1 = vsample
  t2 = ssample
  t3 = betasample
  t4 = tausample
  t5 = eta2sample
  t6 = pisample
  t7 = SS
  dat = list(t1=t1,t2=t2,t3=t3,t4=t4,t5=t5,t6=t6,t7=t7)
  return(dat)
}
