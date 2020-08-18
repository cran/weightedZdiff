\name{zdifference_var}
\alias{zdifference_var}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
z-difference for variance of continuous variable
}
\description{
The function calculates the weighted z-Difference for a continuous variable (x) with binaryreference variable (ref) a}
\usage{
zdifference_var(x,ref,w=NULL,na.rm=TRUE,r)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
The continuous variable for which the weighted z-Difference should be calculated.
}
  \item{ref}{
The binary reference variable as a vector.
}
  \item{w}{
The weights to calculate the weighted binary z-Difference
}
\item{na.rm}{
Should NAs be removed or not. If NAs exists in dataset and na.rm=FALSE then an error will occure.}
\item{r}{
digits to round the returned value, default is 2}
}
\value{
The function returns the calculated z-Difference as a numeric value.
}
\author{
Tim Filla}
%% ~Make other sections like Warning with \section{Warning }{....} ~
\examples{
#generate the data. The weights are taken from uniform distribution and the
#values of x are generated from a bernoulli distribution with success rate 0.3.
#The reference variable is chosen from a bernoulli distribution with success rate 0.8.
ref<-sample(0:1,1000,replace=TRUE,prob=c(0.2,0.8))
erg<-unlist(lapply(1:1000,function(z){
  w<-runif(1000)
  x<-rnorm(1000,1,0.3)
  zdifference_var(x,ref,w)
}))
hist(erg,breaks=50,main="z-difference for continuous data")
plot(seq(0.005,0.97,0.01),quantile(erg,seq(0.005,0.97,0.01)),type="l",
lwd=3,xlab=c("quantile"),ylab=c("x-value"))
points(seq(0.005,0.97,0.01),qnorm(seq(0.005,0.97,0.01)),col="red",type="l",lwd=3,lty=2)
legend("topleft",legend=c("N(0,1) distribution","sample distribution"),lty=c(2,1),
lwd=c(3,3),col=c("red","black"),cex=1.3)
}
