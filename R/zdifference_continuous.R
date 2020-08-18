zdifference_continuous<-function(x,ref,w=NULL,na.rm=TRUE,r=2){
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
  if(class(x)=="integer"){
    x<-as.numeric(as.character(x))
  }
  if(class(x)=="factor"){
    x<-as.numeric(x)
  }
  x0<-x[which(ref==sort(unique(ref))[1])]
  x1<-x[-which(ref==sort(unique(ref))[1])]
  w0<-w[which(ref==sort(unique(ref))[1])]/sum(w[which(ref==sort(unique(ref))[1])])
  w1<-w[which(ref==sort(unique(ref))[2])]/sum(w[which(ref==sort(unique(ref))[2])])
  xw0<-sum(w0*x0,na.rm=na.rm)/sum(w0,na.rm=na.rm)
  xw1<-sum(w1*x1,na.rm=na.rm)/sum(w1,na.rm=na.rm)
  var0<-sum(w0^2,na.rm=na.rm)*var(x0,na.rm=na.rm)
  var1<-sum(w1^2,na.rm=na.rm)*var(x1,na.rm=na.rm)
  erg<-(xw0 - xw1)/sqrt(var0+var1)
  return(round(erg,r))
}
