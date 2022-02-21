center.data<-function(data, predictors){
  
     center.data.output<- sapply(predictors, function(x){
      
     center.data<-data[,colnames(data) %in% x]
  
     center.data.output<-(center.data-mean(center.data, na.rm=TRUE))/
     sd(center.data, na.rm=TRUE)
      
      return(center.data.output)
    
  })
    new.data<-as.data.frame(cbind(data[,!colnames(data) %in% predictors]),
                            center.data.output)
  return(new.data)
}