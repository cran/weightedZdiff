\name{zdifference_nominal}
\alias{zdifference_nominal}
\title{
z-difference for nominal variables.
}
\description{
The function calculates the nominal weighted z-Difference for a binary    reference variable (ref) and a nominal variable (x)}
\usage{
zdifference_nominal(x,ref,w=NULL,na.rm=TRUE,norma=TRUE,r=2)
}
\arguments{
  \item{x}{
The continuous variable for which the weighted z-Difference should be calculated.
}
  \item{ref}{
The binary reference variable as a vector.
}
  \item{w}{
The weights to calculate the weighted continuous z-Difference
}
  \item{na.rm}{
Should NAs be removed or not. If NAs exists in dataset and na.rm=FALSE then an error will occure.}
 \item{norma}{
 If norma = TRUE the weighted z-Difference has a standard Gaussian distribution. If norma = FALSE the resulting distribution is chi squared with #status -1 as degree of freedom.
 }
  \item{r}{
digits to round the returned value, default is 2
}
}
\value{
The function returns the calculated z-Difference as a numeric value.
}
\author{
Tim Filla}
%% ~Make other sections like Warning with \section{Warning }{....} ~
\examples{
#generate data. The weights are taken from uniform distribution and the
#values of x are generated from a multinomial distribution with success
#rate (0.2,0.2,0.3,0.15,0.15) for the five different status.The reference
#variable is chosen from a bernoulli distribution with success rate 0.8.
ref<-sample(1:0,1000,replace=TRUE,prob=c(0.2,0.8))
erg<-unlist(lapply(1:1000,function(z){
  w<-runif(1000)
  x<-sample(0:4,1000,replace=TRUE,prob=c(0.2,0.2,0.3,0.15,0.15))
  zdifference_nominal(x,ref,w,norma=TRUE)
}))
hist(erg,breaks=50,main="z-difference for nominal data")
plot(seq(0.005,0.97,0.01),quantile(erg,seq(0.005,0.97,0.01)),type="l",lwd=3)
points(seq(0.005,0.97,0.01),qnorm(seq(0.005,0.97,0.01)),col="red",type="l",lwd=2,lty="dashed")
}
