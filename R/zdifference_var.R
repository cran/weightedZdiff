zdifference_var<-function(x,ref,w=NULL,na.rm=TRUE,r=2){
  if(is.null(w)){
    w<-rep(1,length(x))
  }
  if(na.rm){
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
  x0<-scale(x[which(ref== c(sort(unique(ref))[1]))],scale=FALSE)
  x1<-scale(x[which(ref==c(sort(unique(ref))[2]))],scale=FALSE)
  w0<-w[which(ref==sort(unique(ref))[1] & is.na(x)==FALSE)]/sum(w[which(ref==sort(unique(ref))[1])])
  w1<-w[which(ref==sort(unique(ref))[2]& is.na(x)==FALSE)]/sum(w[which(ref==sort(unique(ref))[2])])
  resultat<-function(w,x){
    return(-3*((sum(w^4*sum(w))/(sum(w)^4)*mean(x^4))+(4*mean(x^3)*mean(x)*(sum(w^3-w^4)))+
                 (3*(var(x)+mean(x)^2)^2*(sum(w^2*sum(w^2))-sum(w^4)))+(6*(var(x)+mean(x)^2)*mean(x)^2*(sum(w^4-w^3)+sum(w^2*sum(-w^2+w))+sum(w^4-w^3)))+
                 (mean(x)^4*(sum(w*(1-w)^2*(1-w))+sum(w*sum(w^3-2*w^2))-sum(w^4-2*w^3)+sum(2*w^2*sum(w^2))-sum(2*w^4)+sum(w*sum(w*sum(-w^2)))-sum(-w^4)-sum(w^2*sum(-w^2))+2*sum(-w^4)-sum(w*sum(-w^3))-sum(-w^3*(1-w)))))+2*(
                   (sum(w^3*sum(w))*mean(x^4)/sum(w)^2)+(2*mean(x^3)*mean(x)*sum(w^2-w^3))+
                     ((var(x)+mean(x)^2)^2*(sum(w*sum(w^2))-sum(w^3)))+((var(x)+mean(x)^2)*mean(x)^2*(sum(w*sum(w-w^2))+sum(w^3-w^2)+sum(w^3 -w^2))))-4*(
                       (sum(w^3)/sum(w)*mean(x^4))+2*(-sum(w^3)/sum(w)*mean(x^3)*mean(x)+sum(w^2*sum(w))*mean(x^3)*mean(x))+
                         (-sum(w^3)/sum(w)*(mean(x)^2+var(x))^2+sum(w^2*sum(w))*(mean(x)^2+var(x))^2)+(var(x)+mean(x)^2)*mean(x)^2*(sum(w*sum(w-w^2))+sum(-w^2+w^3)+sum(w^3-w^2)))+
             (sum(w*sum(w))*(var(x)+mean(x)^2)^2-sum(w^2)*(var(x)+mean(x)^2)^2+sum(w^2)*mean(x^4))+
             4*((sum(w^4)/sum(w)^2*mean(x^4))+(6/sum(w)^2*(sum(w^2*sum(w^2))*(var(x)+mean(x)^2)^2-sum(w^4)*(var(x)+mean(x)^2)^2))+
                  (4/sum(w)^2*(sum(w^3*sum(w))*mean(x)^3*mean(x)-sum(w^4)*mean(x)^3*mean(x))))-
             ((sum(w)*(mean(x)^2+var(x)))+1/2*(-2*sum(w*(1-w))*mean(x)^2-2*sum(w^2)*(mean(x)^2+var(x))))^2)

  }
  round((sum(w0*(x0-sum(w0*x0))^2)/(1-sum(w0^2))-sum(w1*(x1-sum(w1*x1))^2)/(1-sum(w1^2)))/(sqrt(resultat(w0,x0)/(1-sum(w0^2))^2+resultat(w1,x1)/(1-sum(w1^2))^2)),r)
}
