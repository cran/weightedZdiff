zdifference_ordinal<-function(x,ref,w=NULL,na.rm=TRUE,r=10){
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
  x0<-x[which(ref==(sort(unique(ref)))[1])]
  x1<-x[-which(ref==(sort(unique(ref)))[1])]
  w0<-w[which(ref==(sort(unique(ref)))[1])]/sum(w[which(ref==sort(unique(ref))[1])],na.rm=na.rm)
  w1<-w[which(ref==sort((unique(ref)))[2])]/sum(w[which(ref==sort(unique(ref))[2])],na.rm=na.rm)
  x<-c(x0,x1)
  w<-c(w0,w1)
  d<-length(which(ref==sort(unique(ref))[1]))
  ran<-rank(x,ties.method="average")
  t<-unlist(lapply(1:length(unique(ran)),function(i) length(which(ran==unique(ran)[i]))))
  if(length(unique(ran))==length(x)){
    nenner<-sum(w0^2+w1^2)*((length(x)^2-1)/12+(length(x)+1)/12)
  }
  else{
    rs<-unique(ran)
    N<-length(ran)
    cvv<-sum(rs^2*t*(t-1))/(N*(N-1))+(sum(rs*t*sum(rs*t)-rs^2*t^2))/(N*(N-1))-((N+1)/2)^2
    v<-c(-sum(t^3-t)/12/length(ran)+(length(ran)^2-1)/12)
    nenner<-sum(w0^2)*v+(-sum(w0^2))*cvv+sum(w1^2)*v +(-sum(w1^2))*cvv
  }
  return(round((sum(w0*ran[1:d])-sum(w1*ran[-c(1:d)]))/sqrt(nenner),r))
}
