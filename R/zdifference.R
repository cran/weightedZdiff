zdifference<-function(dataset,ref,weights=NULL,standard_weights=FALSE,na.rm=TRUE,binary_variable=NULL,ordinal_variable=NULL
                      ,continuous_variable=NULL,nominal_variable=NULL,r=2,var.est=FALSE,coefvar.est=FALSE,grad=1){
  if(is.null(weights)){
    w<-as.matrix(rep(1,dim(dataset)[1]))
    indw<-1
  }
  else{
    if(any(weights%in%names(dataset))){
      if(any(!weights%in%names(dataset))){
        warning(paste("weight names '",weights[which(!weights%in%names(dataset))],"'not found in the names of the dataset.",sep=""))
      }
      indw<-which(weights%in%names(dataset))
      w<-as.matrix(dataset[,weights[which(weights %in% names(dataset))]])
    }
    else{
      stop(paste("the name ",w," of the reference variable does not exist in the names of the dataset"))
    }
  }
  if(standard_weights==TRUE & !is.null(weights)){
    w<-cbind(w,rep(1,length(w[,1])))
  }
  if(length(unique(w))!=1 & coefvar.est==TRUE){
    coefvar.est=FALSE
    warning("weighted zdifference for coefficient of variation is not available for non unique weights")
  }
  if(length(unique(na.omit(dataset[,ref])))!= 2){
    stop("reference variable is not binary")
  }
  if(any(w<0)){
    stop("negative weights are not allowed")
  }
  reference<-dataset[,ref]
  dataset<-dataset[,-c(which(colnames(dataset)==ref),which(colnames(dataset)%in%weights==TRUE))]
  res<-rep(NA,ncol(dataset))
  if(!is.null(continuous_variable)){
    if(any(!continuous_variable %in% colnames(dataset))){
      stop("At least one of the variables in continuous_variable does not specify a column of data.")
    }
    res[which(colnames(dataset) %in% continuous_variable)]<-"con"
  }
  if(!is.null(binary_variable)){
    if(any(!binary_variable %in% colnames(dataset))){
      stop("At least one of the variables in binary_variable does not specify a column of data.")
    }
    res[which(colnames(dataset) %in% binary_variable)]<-"bin"
  }
  if(!is.null(ordinal_variable)){
    if(any(!ordinal_variable %in% colnames(dataset))){
      stop("At least one of the variables in ordinal_variable does not specify a column of data.")
    }
    res[which(colnames(dataset) %in% ordinal_variable)]<-"ordi"
  }

  if(!is.null(nominal_variable)){
    if(any(!nominal_variable %in% colnames(dataset))){
      stop("At least one of the variables in nominal_variable does not specify a column of data.")
    }
    res[which(colnames(dataset) %in% nominal_variable)]<-"nom"
  }
  if(any(is.na(res))){
    for(i in 1:ncol(dataset)){
      if(is.na(res[i])){
        if(length(unique(dataset[,i]))==2){
          res[i]<-"bin"
        }
        else if(length(unique(dataset[,i]))< 10 | class(dataset[,i])=="factor"){
          res[i]<-"nom"
        }
        else if(length(unique(dataset[,i])) >= 10 ){
          res[i]<-"con"
        }
      }
    }
    print("User did not specify all variable types. Please check if the variables are correctly assigned and if not change the second column.
          'con' means continuous, 'bin' binary, 'ordi' ordinal and 'nom' nominal variables.")
    zdiffoverview<-data.frame("names"=colnames(dataset),"type"=res)
    res<-edit(zdiffoverview)[,2]
    while(any(!res%in%c("nom","con","bin","ordi"))){
      print("User did not specify all variable types. Please check if the variables are correctly assigned and if not changethe second column.
            'con' means continuous, 'bin' binary, 'ordi' ordinal and 'nom' nominal variables.")
      res<-edit(zdiffoverview)[,2]
    }
  }
  zdiff<-matrix(NA,dim(dataset)[2],dim(w)[2])
  for(i in 1:dim(w)[2]){
    if(length(unique(w[,i]))==1){
      zdiff[,i]<-unlist(lapply(1:dim(dataset)[2],function(z){
        if(res[z]=="con"){
          return(zdifference_continuous(dataset[,z],reference,w=w[,i],r=r))
        }
        else if(res[z]=="bin"){
          return(zdifference_binary(dataset[,z],reference,w=w[,i],r=r))
        }
        else if(res[z]=="ordi"){
          return(zdifference_ordinal(dataset[,z],reference,r=r))
        }
        else if(res[z]=="nom"){
          return(zdifference_nominal(dataset[,z],reference,w=w[,i],r=r))
        }
      }))
    }
    else{
      zdiff[,i]<-unlist(lapply(1:dim(dataset)[2],function(z){
        if(res[z]=="con"){
          return(zdifference_continuous(dataset[,z],reference,w[,i],r=r))
        }
        else if(res[z]=="bin"){
          return(zdifference_binary(dataset[,z],reference,w[,i],r=r))
        }
        else if(res[z]=="ordi"){
          return(zdifference_ordinal(dataset[,z],reference,w[,i],r=r))
        }
        else if(res[z]=="nom"){
          if(any(table(dataset[,z])/length(dataset[,z])<0.06)){
            summe<-0
            for(i in 1:length(na.omit(unique(dataset[,z])))){
              newvar<-ifelse(dataset[,z]==unique(dataset[,z])[i],1,0)
              summe<-summe+zdifference_binary(newvar,reference,w[,i],r=r)
            }
            return(summe)
          }
          else{
            return(zdifference_nominal(dataset[,z],reference,w[,i],r=r))
          }
        }
      }))
    }
  }
  if(grad!=1){
    zdiffgrad<-list()
    for(i in 1:dim(w)[2]){
      zdiffgrad[[i]]<-matrix(NA,length(which(res=="con")),grad-1)
      for(g in 2:grad){
        if(length(unique(w[,i]))==1){
          zdiffgrad[[i]][,g-1]<-unlist(lapply(1:dim(dataset)[2],function(z){
            if(res[z]=="con"){
              return(zdifference_continuous(dataset[,z]^g,reference,w=w[,i],r=r))
            }
          }))
        }
        else{
          zdiffgrad[[i]][,g-1]<-unlist(lapply(1:dim(dataset)[2],function(z){
            if(res[z]=="con"){
              return(zdifference_continuous(dataset[,z]^g,reference,w=w[,i],r=r))
            }
          }))
        }
      }
    }
  }
  if(var.est){
    zdiffvar<-matrix(NA,length(which(res=="con")),dim(w)[2])
    for(i in 1:dim(w)[2]){
      zdiffvar[,i]<-unlist(lapply(1:dim(dataset)[2],function(z){
        if(res[z]=="con"){
          return(zdifference_var(dataset[,z],reference,w[,i],r=r))
        }
      }))
    }
  }
  if(coefvar.est){
    zdiffcoefvar<-matrix(NA,length(which(res=="con")),dim(w)[2])
    for(i in 1:dim(w)[2]){
      zdiffcoefvar[,i]<-unlist(lapply(1:dim(dataset)[2],function(z){
        if(res[z]=="con"){
          return(zdifference_coefvar(dataset[,z],reference,w[,i],r=r))
        }
      }))
    }
  }
  if(grad!=1){
    gradsumme<-unlist(lapply(1:dim(dataset)[2],function(z){
      if(res[z]=="con"){
        summe<-0
        for(i in 2:grad){
          summe<-summe+abs(zdifference_continuous(dataset[,z]^i,reference,w,r=r))
        }
        return(summe)
      }
    }))
  }
  else{
    gradsumme<-0
  }
  if(grad==1){
    moment=c(rep(1,length(res)))
    moment[which(res!="con")]<-""
    zdiffmatrix<-matrix(NA,length(c(zdiff[,1])),dim(w)[2])
    for(i in 1:dim(w)[2]){
      zdiffmatrix[,i]<-zdiff[,i]
    }
  }
  if(grad!=1){
    moment=c(rep(1,length(res)),rep(2:grad,each=length(which(res=="con"))))
    moment[which(res!="con")]<-""
    print(zdiffgrad)
    print(zdiff)
    zdiffmatrix<-matrix(NA,length(zdiffgrad[[1]])+length(zdiff)/dim(w)[2],dim(w)[2])
    print(zdiffmatrix)
    for(i in 1:dim(w)[2]){
      zdiffmatrix[,i]<-cbind(c(zdiff[,i],zdiffgrad[[i]]))
    }
  }

  res<-as.character(res)
  variable<-data.frame("variable"=c(rep(colnames(dataset),1),rep(colnames(dataset)[which(res=="con")],grad-1)),"zdifference"=zdiffmatrix,"moment"=moment)
  if(standard_weights==TRUE & !is.null(weights)){
    colnames(variable)<-c("variable",c(weights[indw],"standard weights"),"moment")
  }
  else if(is.null(weights)){
    colnames(variable)<-c("variable","unweighted","moment")
  }
  else{
    colnames(variable)<-c("variable",weights[indw],"moment")
  }

  if(var.est==TRUE & coefvar.est==TRUE){
    var.variable<-data.frame("variable"=colnames(dataset)[which(res=="con")],
                             "zdifference"=zdiffvar,"zdiffcoefvar"=zdiffcoefvar,"treated"=rep("con",length(which(res=="con"))))
    return(list("result"= variable,"variance_result"=var.variable, "squaredzdiff"=colSums(abs(zdiffmatrix)^2)))
  }
  else if(var.est==TRUE & coefvar.est==FALSE){
    var.variable<-data.frame("variable"=colnames(dataset)[which(res=="con")],
                             "zdifference"=zdiffvar,"treated"=rep("con",length(which(res=="con"))))
    return(list("result"= variable,"variance_result"=var.variable, "squaredzdiff"=colSums(abs(zdiffmatrix)^2)))
  }
  else if(var.est==FALSE & coefvar.est==TRUE){
    var.variable<-data.frame("variable"=colnames(dataset)[which(res=="con")],
                             "zdiffcoefvar"=zdiffcoefvar,"treated"=rep("con",length(which(res=="con"))))
    return(list("result"= variable,"variance_result"=var.variable, "squaredzdiff"=colSums(abs(zdiffmatrix)^2)))
  }
  else{
    return(list("result"= variable, "squaredzdiff"=colSums(abs(zdiffmatrix)^2)))
  }
  }
