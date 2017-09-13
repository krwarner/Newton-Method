
#########################
# Load in data set:
#########################

setwd("~/Desktop/Data")
data <- read.table("Example71.txt")

#########################
# Plot of raw data
#########################
plot(data$x,data$y,xlab="Tumor size (cm), X", ylab="Lymph node metastasis, Y",main="Raw Data",ylim=c(-.2,1.2),xlim=c(0,12))
abline(h=0,lty=3)
abline(h=1,lty=3)
abline(v=min(data$x),lty=2,col=1)
abline(v=max(data$x),lty=2,col=1)

#########################
# Estimate logistic model usinf glm()
#########################

model <- glm(y~x,data=data,family=binomial(link = "logit"))
linear.pred <- predict(model,newdata=data.frame(x=7))
probs <- exp(linear.pred)/(1+exp(linear.pred))
probs


#########################
# Plot of raw data and estimated model
#########################

plot(data$x,data$y,xlab="Tumor size (cm), X", ylab="Lymph node metastasis, Y",main="Raw Data",ylim=c(-.2,1.2),xlim=c(0,12))
abline(h=0,lty=3)
abline(h=1,lty=3)
abline(v=min(data$x),lty=2,col=1)
abline(v=max(data$x),lty=2,col=1)
x.pred <- seq(-1,13,by=.1)
linear.pred <- predict(model,newdata=data.frame(x=x.pred))
probs <- exp(linear.pred)/(1+exp(linear.pred))
lines(x.pred,probs,col="purple")


#########################
# Maximum likelihood estimation 
#########################

#########################
# Log-likelihood function (negative)
#########################

logistic.Neg.LL <- function(b,data=data) {
  
  b0 <- b[1]
  b1 <- b[2]
  x <- data$x
  y <- data$y
  p.i <- exp(b0+b1*x)/(1+exp(b0+b1*x))
  return(-sum(dbinom(y,size=1,prob=p.i,log=TRUE)))
  
}

logistic.Neg.LL(b=c(0,0),data=data)
logistic.Neg.LL(b=c(-3,.6),data=data)

mode(x)


#########################
# grad.descent from class
#########################

library(numDeriv)

Newton.Method <- function(f, x0, max.iter = 200, stopping.deriv = 0.01, ...) {
  
  n    <- length(x0)
  xmat <- matrix(0, nrow = n, ncol = max.iter)
  xmat[,1] <- x0
  
  for (k in 2:max.iter) {
    # Calculate the gradient
    grad.cur <- grad(f, xmat[ ,k-1], ...) 
    
    #Calc Hessian
    hess.cur <- hessian(f, xmat[, k-1], ...)
    
    # Should we stop?
    if (all(abs(grad.cur) < stopping.deriv)) {
      k <- k-1; break
    }
    
    # Move in the opposite direction of the grad
    xmat[ ,k] <- xmat[ ,k-1] - solve(hess.cur) %*% grad.cur
  }
  
  xmat <- xmat[ ,1:k] # Trim
  return(list(x = xmat[,k],
              xmat = xmat, 
              k = k, 
              minimum=f(xmat[,k],...)
  )
  )
}

#########################
# Optimization using grad.descent()
#########################

x0 <- c(0,0)
nm <- Newton.Method(logistic.Neg.LL,x0,max.iter = 2000, stopping.deriv = 0.001,data=data)
nm$minimum
nm$k
nm$x

#########################
# Optimization using nlm()
#########################


nml.opt <- nlm(logistic.Neg.LL,c(0,0),data=data)
nml.opt




