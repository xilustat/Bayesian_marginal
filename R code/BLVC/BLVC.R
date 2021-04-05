BLVC = function(e,c,x,w, y, max.steps) 
{
  
  n <- length(y)
  q1 <- ncol(c)
  q2 <- ncol(e)
  a1 = 1 #lambda2_c (prior)~ gamma(a1,b1)
  b1 = 1
  a2 = 1 #lambda2_e (prior)~ gamma(a2,b2)
  b2 = 1
  
  beta = 1 # coeffecient of gene
  b = rep(1,q1) # coeffecient of clinical factor
  alpha = rep(1,q2) # coeffecient of environmental factor
  eta = rep(1,q2) # coeffecient of GXE interaction 
  invtau2_c = 1
  invtau2_e = rep(1,q2)
  lambda2_c = 1
  lambda2_e = 1
  sigma2 = 1 
  
  betaSamples <- matrix(0,max.steps,1)
  bSamples <- matrix(0,max.steps,q1)
  alphaSamples <- matrix(0,max.steps,q2)
  etaSamples <- matrix(0, max.steps, q2)
  sigma2Samples <- matrix(0, max.steps,1)
  invTau2_cSamples <- matrix(0, max.steps, 1)
  invTau2_eSamples <- matrix(0, max.steps, q2)
  lambda2_cSamples = matrix(0, max.steps,1)
  lambda2_eSamples <- matrix(0, max.steps,1)
  
  for(k in 1:max.steps)
  {
  
    #sample b
    invsig1 = solve(diag(1,nrow=q1))
    B = invsig1+t(c)%*%c/sigma2
    varcov1 = solve(B)
    res1 = y-e%*%alpha-x*beta-w%*%eta
    mean1 = varcov1%*%t(t(res1)%*%c/sigma2)
    b = mvrnorm(1,mean1,varcov1)
    
    bSamples[k,] <- b
    
    #sample alpha
    invsig2 = solve(diag(1,nrow=q2))
    B1 = invsig2+t(e)%*%e/sigma2
    varcov2 = solve(B1)
    res2 = y-c%*%b-x*beta-w%*%eta
    mean2 = varcov2%*%t(t(res2)%*%e/sigma2)
    alpha = mvrnorm(1,mean2,varcov2)
    
    alphaSamples[k,] <- alpha
    
    # sample beta
    sigB = 1/(t(x)%*%x+invtau2_c)
    varB = sigma2*sigB
    meanB = sigB*t(x)%*%(y-c%*%b-e%*%alpha-w%*%eta)
    beta = rnorm(1,meanB,varB)
    
    betaSamples[k,] <- beta
    
    # sample eta
    sigeta = solve(t(w)%*%w + diag(invtau2_e, nrow = q2))
    varcov3 = sigma2*sigeta
    mean3 = sigeta%*%t(w)%*%(y-c%*%b-e%*%alpha-x*beta)
    eta = mvrnorm(1, mean3, varcov3)
    
    etaSamples[k,] <- eta
    
    # sample sigma2
    shape = (n+1+q2)/2
    residue = y-c%*%b-e%*%alpha-x*beta-w%*%eta
    scale = t(residue)%*%residue+beta^2*invtau2_c+t(eta)%*%diag(invtau2_e, nrow = q2)%*%eta
    sigma2 = 1/rgamma(1, shape, rate = scale)
    
    sigma2Samples[k,] <- sigma2
    
    #sample invtau2_c
    muPrime1 = sqrt(lambda2_c * sigma2 / beta^2)
    lambdaPrime1 = lambda2_c
    invtau2_c = rinvgauss(1, m=muPrime1, s=lambdaPrime1)
    
    invTau2_cSamples[k, ] <- invtau2_c
    
    #sample invtau2_e
    muPrime2 <- sqrt(lambda2_e * sigma2 / eta^2)
    lambdaPrime2 <- lambda2_e
    #invtau2_e <- rep(0, q2)
    for (i in 1:q2) {
      invtau2_e[i] <- rinvgauss(1, m=muPrime2[i], s=lambdaPrime2)
    }
    
    invTau2_eSamples[k, ] <- invtau2_e
    
    # sample lambda2_c
    shape1 = a1 + 1
    rate1 = b1 + (1/invtau2_c)/2
    lambda2_c <- rgamma(1, shape1, rate = rate1)
    
    lambda2_cSamples[k,] <- lambda2_c
    
    # sample lambda2_e
    shape2 = a2 + q2
    rate2 = b2 + sum(1/invtau2_e)/2
    lambda2_e <- rgamma(1, shape2, rate = rate2)
    
    lambda2_eSamples[k,] <- lambda2_e
    
    
  }
  
  t1 = betaSamples
  t2 = etaSamples
  t3 = sigma2Samples
  t4 = lambda2_cSamples
  t5 = lambda2_eSamples
  t6 = invTau2_cSamples
  t7 = invTau2_eSamples
  t8 = bSamples
  t9 = alphaSamples
    
  
  dat = list(t1=t1,t2=t2,t3=t3,t4=t4,t5=t5,t6=t6,t7=t7,t8=t8,t9=t9)
  return(dat)
  
}

