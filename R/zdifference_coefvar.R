zdifference_coefvar<-function(x,ref,na.rm=TRUE,r=2){
  if(length(unique(na.omit(ref)))!= 2){
    stop("reference variable is not binomial")
  }
  if(class(x)=="integer"){
    x<-as.numeric(as.character(x))
  }
  if(class(x)=="factor"){
    x<-as.numeric(x)
  }
  if(na.rm){
    exclna<-which(is.na(x)| is.na(ref))
    if(length(exclna)!=0){
      x<-x[-exclna]
      ref<-ref[-exclna]
    }
  }
  else{
    if(any(is.na(x)| is.na(ref)))
      warning("NA in the data and NAs not removed")
  }
  x0<-x[which(ref==sort(unique(ref))[1])]
  x1<-x[-which(ref==sort(unique(ref))[1])]
  n0<-length(x0)
  n1<-length(x1)
  mean0<-mean(x0)
  mean1<-mean(x1)
  var0<-var(x0)
  var1<-var(x1)
  c0<-sqrt(var0)/mean0
  c1<-sqrt(var1)/mean1
  cv<-((n0-1)*c0+(n1-1)*c1)/(n0+n1-2)
  res<-(c0-c1)/sqrt((1/(n0-1)+1/(n1-1))*cv^2*(0.5+cv^2))
  return(round(res,r))
}
