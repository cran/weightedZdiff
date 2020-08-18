\name{zdifference_ordinal}
\alias{zdifference_ordinal}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
weighted z-difference for ordinal variables
}
\description{
The function calculates the ordinal weighted z-Difference for a binary    reference variable (ref) and an ordinal variable (x)}
\usage{
zdifference_ordinal(x,ref,w=NULL,na.rm=TRUE,r=10)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
The ordinal variable for which the weighted z-Difference should be calculated.
}
  \item{ref}{
The binary reference variable as a vector.
}
  \item{w}{
The weights to calculate the weighted ordinal z-Difference
}
  \item{r}{
digits to round the returned value, default is 2
}
  \item{na.rm}{
%%     ~~Describe \code{na.rm} here~~
}
}
\value{
The function returns the calculated z-Difference as a numeric value.
}
\references{
For standard z-difference (unweighted)
https://pubmed.ncbi.nlm.nih.gov/23972521/}
\author{
Tim Filla}
%% ~Make other sections like Warning with \section{Warning }{....} ~
\examples{
#generate the data. The weights are taken from uniform distribution and the
#variable x has 5 different status with probability of beeing in status i is
#given by:0.1,0.2,0.3,0.3,0.1. The reference variable
#is chosen from a bernoulli distribution with success rate 0.8.
ref<-sample(0:1,1000,replace=TRUE,prob=c(0.2,0.8))
erg<-unlist(lapply(1:1000,function(z){
  w<-runif(1000)
  x<-sample(1:5,1000,replace=TRUE,prob=c(0.1,0.2,0.3,0.3,0.1))
  zdifference_ordinal(x,ref,w)
}))
hist(erg,breaks=50,main="z-difference for continuous data")
plot(seq(0.005,0.97,0.01),quantile(erg,seq(0.005,0.97,0.01)),type="l",
lwd=3,xlab=c("quantile"),ylab=c("x-value"))
points(seq(0.005,0.97,0.01),qnorm(seq(0.005,0.97,0.01)),col="red",type="l",lwd=3,lty=2)
legend("topleft",legend=c("N(0,1) distribution","sample distribution"),lty=c(2,1),
lwd=c(3,3),col=c("red","black"),cex=1.3)
}
