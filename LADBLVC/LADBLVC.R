LADBLVC = function(e,c,x,w, y, max.steps) 
{
  
  n <- length(y)
  q1 <- ncol(c)
  q2 <- ncol(e)
  a = 0.1
  b = 0.1
  c1 = 0.1
  d1 = 0.1
  c2 = 0.1
  d2 = 0.1
  
  v = rep(1,n)
  s1 = 1
  s2 = rep(1,q2)
  beta = 1
  eta = rep(1,q2)
  b0 = 1
  alpha0 = 1
  tau = 1
  theta = 0.5
  e2 = 2/(theta*(1-theta))
  eta2_1 = 1
  eta2_2 = 1
  b = rep(1,q1)
  alpha = rep(1,q2)
  
  betaSamples <- matrix(0,max.steps,1)
  bSamples <- matrix(0,max.steps,q1)
  alphaSamples <- matrix(0,max.steps,q2)
  etaSamples <- matrix(0, max.steps, q2)
  vSamples <- matrix(0,max.steps,n)
  s1Samples <- matrix(0,max.steps,1)
  s2Samples <- matrix(0,max.steps,q2)
  tauSamples = matrix(0,max.steps,1)
  eta2_1Samples = matrix(0,max.steps,1)
  eta2_2Samples = matrix(0,max.steps,1)
  
  for(k in 1:max.steps)
  {
    
    #sample b
    for(j in 1:q1)
      {
      A = c[,j]^2/v
      invsigma2 = tau*sum(A)/e2+1/b0
      sigma2 = 1/invsigma2
      y_j = as.vector(y -c[,-j]%*%b[-j]-e%*%alpha-x*beta-w%*%eta)
      B = y_j*c[,j]/v
      mu = tau*sum(B)*sigma2/e2
      b[j] = rnorm(1,mean=mu,sd=sqrt(sigma2))
    }
    
    bSamples[k,] <- b
    
    # sample alpha
    for(j in 1:q2)
    {
      A = e[,j]^2/v
      invsigma2 = tau*sum(A)/e2+1/alpha0
      sigma2 = 1/invsigma2
      y_j = as.vector(y -e[,-j]%*%alpha[-j]-c%*%b-x*beta-w%*%eta)
      B = y_j*e[,j]/v
      mu = tau*sum(B)*sigma2/e2
      alpha[j] = rnorm(1,mean=mu,sd=sqrt(sigma2))
    }
    
    alphaSamples[k,] <- alpha
    
    #sample v
    for(i in 1:n)
    {
      res = y-e%*%alpha-c%*%b-x*beta-w%*%eta
      lambda = 2*tau
      mu = sqrt(lambda*e2/(tau*res[i]^2))
      flag=1
      while(flag){
        inv_v= rinvGauss(1,nu=mu,lambda=lambda)
        flag = any(inv_v<=0)|any(is.na(inv_v))|any(is.infinite(inv_v))
      }
      v[i]=1/inv_v
    }
    
    vSamples[k,] = v
    
    #sample s1
    lambda_1 = eta2_1
    mu_1 = sqrt(lambda_1/beta^2)
    invs = rinvGauss(1, nu=mu_1, lambda=lambda_1)
    s1=1/invs
    
    s1Samples[k,] = s1
    
    #sample s2
    for(j in 1:q2){
        lambda_2 = eta2_2
        mu_2 = sqrt(lambda_2/eta[j]^2)
        flag=1
        while(flag){
          invs2 = rinvGauss(1,nu=mu_2,lambda=lambda_2)
          flag = (invs2<=0)|(is.na(invs2))|(is.infinite(invs2))
        }
        s2[j]=1/invs2
      }
    
    s2Samples[k,] = s2
    
    #sample beta
    A = x^2/v
    invsigma2 = tau*sum(A)/e2+1/s1
    sigma2 = 1/invsigma2
    y_j = as.vector(y-e%*%alpha-c%*%b-w%*%eta)
    B = y_j*x/v
    mu = tau*sum(B)*sigma2/e2
    beta = rnorm(1,mean=mu,sd=sqrt(sigma2))
   
    betaSamples[k,]=beta
    
    #sample eta
    for(j in 1:q2){
      A = w[,j]^2/v
      invsigma2 = tau*sum(A)/e2+1/s2[j]
      sigma2 = 1/invsigma2
      y_j = as.vector(y - e%*%alpha-c%*%b-x*beta-w[,-j]%*%eta[-j])
      B = y_j*w[,j]/v
      mu = tau*sum(B)*sigma2/e2
      eta[j] = rnorm(1,mean=mu,sd=sqrt(sigma2))
    }
    
    etaSamples[k,]=eta
    
    #sample tau
    res = y-c%*%b-e%*%alpha-x*beta-w%*%eta
    vec = res^2/(2*e2*v)+v
    rate = sum(vec)+b
    shape = a+3*n/2
    tau = rgamma(1,shape=shape,rate=rate)
    tauSamples[k,] = tau
    
    #sample eta2_1
    rate_1 = s1/2+d1
    shape_1 = 1+c1
    eta2_1 = rgamma(1,shape=shape_1,rate=rate_1)
    eta2_1Samples[k,] = eta2_1
    
    #sample eta2_1
    rate_2 = sum(s2/2)+d2
    shape_2 = q2+c2
    eta2_1 = rgamma(1,shape=shape_2,rate=rate_2)
    eta2_1Samples[k,] = eta2_1
    
    
  }
  
  
  t1 = betaSamples
  t2 = etaSamples
  t3 = tauSamples
  t4 = s1Samples
  t5 = s2Samples
  t6 = eta2_1Samples
  t7 = eta2_2Samples
  t8 = bSamples
  t9 = alphaSamples
  t10 = vSamples
 
  
  dat = list(t1=t1,t2=t2,t3=t3,t4=t4,t5=t5,t6=t6,t7=t7,t8=t8,t9=t9,t10=t10)
  return(dat)
  
  
  
}

