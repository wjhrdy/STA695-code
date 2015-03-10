BinRatio <- function(n1, x1, n2, x2, step=0.1, initStep=0.1, level=3.84){

if(x1 >= n1) stop("all success?")



MLE <- (x2*n1)/(x1*n2)   #### recall we are to find CI for p2/p1, not p1/p2.
EPS <- .Machine$double.eps

Theta <- function(theta, n1=n1, x1=x1, n2=n2, x2=x2)
  {
       llr <- function(const, n1, x1, n2, x2, theta) {
              npllik1 <- -2*(x1*log((const*n1)/x1) + (n1-x1)*log(((1-const)*n1)/(n1-x1)))
              npllik2 <- -2*(x2*log((const*theta*n2)/x2) + (n2-x2)*log(((1-const*theta)*n2)/(n2-x2)))
        return(npllik1 + npllik2)
        }
       upBD <- min( 1-EPS, 1/theta - EPS)
       temp <- optimize(f = llr, lower = EPS, upper = upBD, n1 = n1, x1 = x1, n2 = n2, x2 = x2,theta = theta)

    cstar <- temp$minimum
    val <- temp$objective
    pvalue <- 1 - pchisq( val, df=1)
    list(`-2LLR` = val, cstar = cstar, Pval=pvalue)
  }

value <- 0
step1 <- step
Lbeta <- MLE - initStep
for( i in 1:8 ) {
while(value < level) {
Lbeta <- Lbeta - step1
value <- Theta(theta = Lbeta, n1=n1, x1=x1, n2=n2, x2=x2)$'-2LLR'
}
Lbeta <- Lbeta + step1
step1 <- step1/10
value <- Theta( theta =Lbeta, n1=n1, x1=x1, n2=n2, x2=x2)$'-2LLR'
}

value1 <- value
value <- 0

Ubeta <- MLE + initStep
for( i in 1:8 ) {
   while(value < level) {
    Ubeta <- Ubeta + step
    value <- Theta(theta=Ubeta, n1=n1, x1=x1, n2=n2, x2=x2 )$'-2LLR'
   }
   Ubeta <- Ubeta - step
   step <- step/10
   value <- Theta(theta=Ubeta, n1=n1, x1=x1, n2=n2, x2=x2)$'-2LLR'
}
return( list(Low=Lbeta, Up=Ubeta, FstepL=step1, FstepU=step,
Lvalue=value1, Uvalue=value) )
}
