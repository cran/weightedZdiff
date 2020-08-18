zdifference_binary<-function(x,ref,w=NULL,na.rm=TRUE,r=2){
  if(is.null(w)){
    w<-rep(1,length(x))
  }
  if(na.rm==TRUE){
    exclna<-which(is.na(x)==TRUE |is.na(w)==TRUE | is.na(ref)==TRUE)
    if(length(exclna)!=0){
      x<-x[-exclna]
      w<-w[-exclna]
      ref<-ref[-exclna]
    }
  }
  else{
    if(any(is.na(x)==TRUE |is.na(w)==TRUE | is.na(ref)==TRUE))
      warning("NA in the data and NAs not removed")
  }
  if(any(w<0)){
    stop("negative weights are not allowed")
  }
  if(length(unique(na.omit(ref)))!= 2){
    stop("reference variable is not binomial")
  }
  if(length(unique(na.omit(x)))!=2){
    stop("variable is not binomial")
  }
  if(class(x)=="integer"){
    x<-as.numeric(as.character(x))
  }
  if(class(x)=="factor"){
    x<-as.numeric(x)
  }

  lx<-unique(sort(x))
  l<-sort(unique(ref))
  x0<-x[which(ref==l[1])]
  x1<-x[which(ref==l[2])]
  w0<-w[which(ref==l[1])]/sum(w[which(ref==l[1])])
  w1<-w[which(ref==l[2])]/sum(w[which(ref==l[2])])
  n0<-length(na.omit(x0))
  n1<-length(na.omit(x1))
  p0<-length(which(x0==lx[1]))/n0
  p1<-length(which(x1==lx[1]))/n1
  res<-(sum(w0*x0,na.rm=na.rm)/sum(w0,na.rm=na.rm)-sum(w1*x1,na.rm=na.rm)/sum(w1,na.rm=na.rm))/sqrt(p0*(1-p0)*sum(w0^2,na.rm=na.rm)/sum(w0,na.rm=na.rm)^2+p1*(1-p1)*sum(w1^2,na.rm=na.rm)/sum(w1,na.rm=na.rm)^2)
  return(round(res,r))
}
