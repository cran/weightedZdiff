\name{zdifference}
\alias{zdifference}
\title{
zdifference for dataset
}
\description{
The function calculates the zdifferences for each variable in a dataset or each column in a matrix (depends on the format of your data). Furthermore the sum of the squared zdifferences is calculated. The variables are set into classes continuous, binary and nominal automatically by the following algorithm. If the variable has only 2 different values its treated as binary. If the variable has more then 9 observations or the class of the variable is factor its treated as nominal and otherwise continuous. The user can specify the type of every variable by hand. }
\usage{
zdifference(dataset,ref,weights=NULL,standard_weights=FALSE,na.rm=TRUE,
binary_variable=NULL,ordinal_variable=NULL,continuous_variable=NULL,nominal_variable=NULL,
r=2,var.est=FALSE,coefvar.est=FALSE,grad=1)
}
\arguments{
  \item{dataset}{
An object of class data.frame or matrix, which contains the variables for which the zDifferences should be calculated and the reference variable in columns.
}
  \item{ref}{
The name of the reference variable, name must be in datasets' names.
}
  \item{weights}{
The name of the variable containing the weights for each observation, name must be in datasets' names.
}
  \item{standard_weights}{
Should the unweighted zdifferences be calculated or not.
}
  \item{na.rm}{
Should NAs be removed or not. If NAs exists in dataset and na.rm=FALSE then an error will occure.}
  \item{binary_variable}{
optional: Names of binomial variables.
}
  \item{ordinal_variable}{
optional: Names of ordinal variables.
}
  \item{continuous_variable}{
optional: Names of continuous variables.
}
  \item{nominal_variable}{
optional: Names of nominal variables.
}
  \item{r}{
Number of digits to round the result.
}
  \item{var.est}{
Should the weighted z-Difference for the variances of continuous variables be reported (TRUE) or not (FALSE)}
  \item{coefvar.est}{
Should the coefficient of variation for continuous variables be reported (TRUE) or not (FALSE)}
  \item{grad}{
The Moments for which to calculate the weighted z-Difference for continuous variables. grad=2 means the first and second moments are calculated.
}
}
\references{
For standard z-difference (unweighted)
https://pubmed.ncbi.nlm.nih.gov/23972521/}
\author{
Tim Filla
}
%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
data(testdata)
#new dataset
zdifference(testdata,"treatment",grad=2,continuous_variable=c("age","meanbp1"),
binary_variable=c("CHF","Cirr","colcan","Coma","lungcan","MOSF","sepsis","female","ARF"))
#generate iptw weights
p<-glm(treatment~.,data=testdata,family="binomial")$fitted.values
testdata$weights<-ifelse(testdata$treatment==0,1/(1-p),1/p)
zdifference(testdata,"treatment",weights="weights",grad=2,
continuous_variable=c("age","meanbp1"),binary_variable=c("CHF","Cirr",
"colcan","Coma","lungcan","MOSF","sepsis","female","ARF"),standard_weights=TRUE)}
