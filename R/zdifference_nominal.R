zdifference_nominal<-function(x,ref,w=NULL,na.rm=TRUE,norma=TRUE,r=2){
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
  if(any(w==0)){
    x<-x[-which(w==0)]
    ref<-ref[-which(w==0)]
    w<-w[-which(w==0)]
  }
  if(length(unique(na.omit(ref)))!= 2){
    stop("reference variable is not binomial")
  }
  x0<-na.omit(x[which(ref==sort(unique(ref))[1])])
  x1<-na.omit(x[-which(ref==sort(unique(ref))[1])])
  w0<-w[which(ref==sort(unique(ref))[1]&is.na(x)==FALSE)]
  w1<-w[which(ref==sort(unique(ref))[2]&is.na(x)==FALSE)]
  l0<-sum(w0)
  l1<-sum(w1)
  x00<-unlist(lapply(1:length(unique(na.omit(x))),function(z){
    sum(w0[which(x0==sort(unique(na.omit(x)))[z])])
  }))
  s00<-unlist(lapply(1:length(unique(na.omit(x))),function(z){
    sum(w0[which(x0==sort(unique(na.omit(x)))[z])]^2)
  }))
  x11<-unlist(lapply(1:length(unique(na.omit(x))),function(z){
    sum(w1[which(x1==sort(unique(na.omit(x)))[z])])
  }))
  s11<-unlist(lapply(1:length(unique(na.omit(x))),function(z){
    sum(w1[which(x1==sort(unique(na.omit(x)))[z])]^2)
  }))
  chi<-sum((l0*x11-l1*x00)^2/(l0^2*s11 + l1^2*s00))
  if(norma==TRUE){
    return(round(sample(c(-1,1),1)*qnorm(0.5+(1-pchisq(chi,length(unique(x))-1,lower.tail=FALSE))/2),r))
  }

  else{
    return(round(chi,r))
  }
}
