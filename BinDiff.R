##########################################################
# The BinDiff function is to created to construct 95% CI #
# for the diffrence between two binomial proportions     #
# using the likelihood ratio test for testing p1-p2      #
##########################################################

BinDiff <- function(n1, x1, n2, x2, step=0.1, initStep=0.1, level=3.84){
  
  # p1 & p2 must be between (0,1)
  if((x1 >= n1) | (x2 >= n2)) stop("all success?")
  if((x1 <= 0) | (x2 <= 0)) stop("all failure?")
  
  
  MLE <- (x1/n1)-(x2/n2)   #### recall we are to find CI for p1-p2, not p2-p1.
  EPS <- .Machine$double.eps
  
  Theta <- function(theta, n1=n1, x1=x1, n2=n2, x2=x2)
  {
    llr <- function(const, n1, x1, n2, x2, theta) {
      # Let theta=p1-p2 ..... (1)
      # Liklehood ratio test for the first sample considering p1=const
      npllik1 <- -2*(x1*log((const*n1)/x1) + (n1-x1)*log(((1-const)*n1)/(n1-x1))) 
      # Liklehood ratio test for the second sample considering (from (1) above) p2=const-theta, and p1=const
      npllik2 <- -2*(x2*log(((const-theta)*n2)/x2) + (n2-x2)*log(((1-(const-theta))*n2)/(n2-x2)))
      # Recall that -2*(the log of the likelihoods ration)~ Chi squared with df=1
      return(npllik1 + npllik2)
    }
    # We have (0 < p1, p2 < 1) & hence (-1 < theta < 1) 
    # Here we are optimizing over p1
    upBD <- min( 1-EPS, 1+theta - EPS) # upper bound for p1 (const)
    loBD <- max( EPS , theta+EPS) # lower bound for p2 (const-theta)
    temp <- optimize(f = llr, lower = loBD, upper = upBD, n1 = n1, x1 = x1, n2 = n2, x2 = x2,theta = theta)
    
    cstar <- temp$minimum
    val <- temp$objective
    pvalue <- 1 - pchisq( val, df=1)
    list(`-2LLR` = val, cstar = cstar, Pval=pvalue)
  }
  
  value <- 0
  step1 <- step
  Lbeta <- MLE
  Ubeta <- MLE
  for( i in 1:8 ) {
    Lbeta <- Lbeta - step1
    while((value < level) & ( Lbeta > -1)) {
      value <- Theta(theta = Lbeta, n1=n1, x1=x1, n2=n2, x2=x2)$'-2LLR'
      Lbeta <- Lbeta - step1
    }
    Lbeta <- Lbeta + step1
    step1 <- step1/10
    if (Lbeta != MLE) {
      value <- Theta( theta =Lbeta, n1=n1, x1=x1, n2=n2, x2=x2)$'-2LLR'
    }
    
  }
  
  value1 <- value
  value <- 0
  
  
  for( i in 1:8 ) {
    Ubeta <- Ubeta + step
    while((value < level) & (Ubeta < 1)) {
      value <- Theta(theta=Ubeta, n1=n1, x1=x1, n2=n2, x2=x2 )$'-2LLR'
      Ubeta <- Ubeta + step
    }
    Ubeta <- Ubeta - step
    step <- step/10
    if (Ubeta != MLE){
      value <- Theta(theta=Ubeta, n1=n1, x1=x1, n2=n2, x2=x2)$'-2LLR'
    }
    
  }
  return( list(Low=Lbeta, Up=Ubeta, FstepL=step1, FstepU=step,
               Lvalue=value1, Uvalue=value) )
}

