BinDiff = function(n1, x1, n2, x2, step=0.1, level=3.84){
  
  if(x1 >= n1) stop("all success?")
  # we are to find CI for p1-p2
  MLE = (x1/n1)-(x2/n2)   
  EPS = .Machine$double.eps
  
  Theta = function(theta, n1=n1, x1=x1, n2=n2, x2=x2)
  {
    llr = function(const, n1, x1, n2, x2, theta) {
      # parameterize p1 as const
      npllik1 = -2*(x1*log((const*n1)/x1) + (n1-x1)*
                      log(((1-const)*n1)/(n1-x1)))
      # parameterize p2 as const-p1
      npllik2 = -2*(x2*log(((const-theta)*n2)/x2) + 
                      (n2-x2)*log(((1-(const-theta))*n2)/(n2-x2))) 
      return(npllik1 + npllik2) ##additive log likelihood
    }
    # we must ensure that 0<const<1 and 0<const-theta<1
    upBD = min(1-EPS, 1+theta-EPS)
    loBD = max(EPS, theta+EPS)
    temp = optimize(f = llr, 
                    lower = loBD, 
                    upper = upBD, 
                    n1 = n1, 
                    x1 = x1, 
                    n2 = n2, 
                    x2 = x2,
                    theta = theta)
    
    cstar = temp$minimum
    val = temp$objective
    pvalue = 1 - pchisq( val, df=1)
    list(`-2LLR` = val, cstar = cstar, Pval=pvalue)
  }
  Lbeta = Ubeta = MLE
  value = 0
  step1 = step
  for( i in 1:8 ) {
    Lbeta = Lbeta - step1
    while((value < level)&(Lbeta > -1)) {
      value = Theta(theta = Lbeta, n1=n1, x1=x1, n2=n2, x2=x2)$'-2LLR'
      Lbeta = Lbeta - step1
    }
    Lbeta = Lbeta + step1
        
    step1 = step1/10
    if (Lbeta!=MLE) {
      value = Theta( theta =Lbeta, n1=n1, x1=x1, n2=n2, x2=x2)$'-2LLR'
    }
  }
  
  value1 = value
  value = 0
  
  for( i in 1:8 ) {
    Ubeta = Ubeta + step
    while((value < level)&(Ubeta < 1)) {
      value = Theta(theta=Ubeta, n1=n1, x1=x1, n2=n2, x2=x2 )$'-2LLR'
      Ubeta = Ubeta + step
    }
    Ubeta = Ubeta - step
    step = step/10
    if (Lbeta!=MLE) {
      value = Theta(theta=Ubeta, n1=n1, x1=x1, n2=n2, x2=x2)$'-2LLR'
    }
  }
  
  return( list(Low=Lbeta, Up=Ubeta, FstepL=step1, FstepU=step,
               Lvalue=value1, Uvalue=value) )
}


BinDiff(100,1,100,99)

