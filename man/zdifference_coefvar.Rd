\name{zdifference_coefvar}
\alias{zdifference_coefvar}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
z-difference for the coefficient of vartiation for normal distributed variables.
}
\description{
The function calculates the coefficient of variation z-Difference for a binary reference variable (ref) and an ordinal variable (x)
}
\usage{
zdifference_coefvar(x,ref,na.rm=TRUE,r=2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
The variable for which the z-Difference should be calculated.
}
  \item{ref}{
The binary reference variable as a vector.
}
  \item{na.rm}{
Should NAs be removed or not. If NAs exists in dataset and na.rm=FALSE then an error will occure.}
  \item{r}{
digits to round the returned value, default is 2
}
}
\value{
The function returns the calculated z-Difference as a numeric value.
}
\references{
https://pubmed.ncbi.nlm.nih.gov/23972521/
}
\author{
Tim Filla}
%% ~Make other sections like Warning with \section{Warning }{....} ~
\examples{
#generate the data.
#variable x has 5 different status with probability of
#beeing in status i is given by:0.1,0.2,0.3,0.3,0.1. #The reference variable
#is chosen from a bernoulli distribution with success #rate 0.8.
ref<-sample(0:1,1000,replace=TRUE,prob=c(0.2,0.8))
erg<-unlist(lapply(1:1000,function(z){
  w<-runif(1000)
  x<-rnorm(1000,25)
  zdifference_coefvar(x,ref)
}))
hist(erg,breaks=50,main="z-difference for continuous data")
plot(seq(0.005,0.97,0.01),quantile(erg,seq(0.005,0.97,0.01)),
type="l",lwd=3,xlab=c("quantile"),ylab=c("x-value"))
points(seq(0.005,0.97,0.01),qnorm(seq(0.005,0.97,0.01)),col="red",type="l",lwd=3,lty=2)
legend("topleft",legend=c("N(0,1) distribution","sample distribution"),lty=c(2,1),
lwd=c(3,3),col=c("red","black"),cex=1.3)
}
