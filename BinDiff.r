BinDiff <- function(n1, x1, n2, x2, step=0.1, initStep=0.1, level=3.84){
  
  if(x1 >= n1) stop("all success?")
  
  MLE <- (x1/n1)-(x2/n2)   #### we are to find CI for p1-p2
  EPS <- .Machine$double.eps
  
  Theta <- function(theta, n1=n1, x1=x1, n2=n2, x2=x2)
  {
    llr <- function(const, n1, x1, n2, x2, theta) {
      npllik1 <- -2*(x1*log((const*n1)/x1) + (n1-x1)*log(((1-const)*n1)/(n1-x1))) ##parameterize p1 as const
      npllik2 <- -2*(x2*log(((const-theta)*n2)/x2) + (n2-x2)*log(((1-(const-theta))*n2)/(n2-x2))) ##parameterize p2 as const-p1
      return(npllik1 + npllik2) ##additive log likelihood
    }
    ##we must ensure that 0<const<1 and 0<const-theta<1
    upBD <- min(1-EPS, 1+theta-EPS)
    loBD <- max(EPS, theta+EPS)
    temp <- optimize(f = llr, lower = loBD, upper = upBD, n1 = n1, x1 = x1, n2 = n2, x2 = x2,theta = theta)
    
    cstar <- temp$minimum
    val <- temp$objective
    pvalue <- 1 - pchisq( val, df=1)
    list(`-2LLR` = val, cstar = cstar, Pval=pvalue)
  }
  
  value <- 0
  step1 <- step
  Lbeta <- MLE - initStep
  if (Lbeta <= -1) stop("try a smaller initial step") ##gives us informative error instead of breaking inside optomize function
  for( i in 1:8 ) {
    while(value < level) {
      Lbeta <- Lbeta - step1
      if (Lbeta <= -1) stop('try a smaller step')
      value <- Theta(theta = Lbeta, n1=n1, x1=x1, n2=n2, x2=x2)$'-2LLR'
    }
    Lbeta <- Lbeta + step1
    step1 <- step1/10
    value <- Theta( theta =Lbeta, n1=n1, x1=x1, n2=n2, x2=x2)$'-2LLR'
  }
  
  value1 <- value
  value <- 0
  
  Ubeta <- MLE + initStep
  if (Ubeta >= 1) stop("try a smaller initial step")
  for( i in 1:8 ) {
    while(value < level) {
      Ubeta <- Ubeta + step
      if (Ubeta >= 1) stop("try a smaller step")
      value <- Theta(theta=Ubeta, n1=n1, x1=x1, n2=n2, x2=x2 )$'-2LLR'
    }
    Ubeta <- Ubeta - step
    step <- step/10
    value <- Theta(theta=Ubeta, n1=n1, x1=x1, n2=n2, x2=x2)$'-2LLR'
  }
  return( list(Low=Lbeta, Up=Ubeta, FstepL=step1, FstepU=step,
               Lvalue=value1, Uvalue=value) )
}


BinDiff(100,1,100,99,initStep=0.01,step=0.001)

