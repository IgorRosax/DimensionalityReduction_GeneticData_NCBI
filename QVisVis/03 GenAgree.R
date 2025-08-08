GenAgree <- function(Config1,Config2,Startk,Endk,DMetric="euclidian",varargin=NULL,Description=NULL){
  #GenAgree
  #Calculates a generalized agreement metric as per France and Carroll (2007)
  #Config 1 - Configuration 1 to compare
  #Config 2 - Configuration 2 to compare
  #Startk - The starting value of k
  #Endk - The end value of k
  #DMetric - The distance metric used to calculate distances
  #DMetricParam - An optional parameter for the distance metric
  #The name of the distance metric for pdist function
  #OutValue - The value of the generalized agreement
  #RBuild - k*4 table giving following values for each k
  # 1 - The value of AR(k)
  # 2 - ER(k) drawn from hypergeometric
  # 3 - Numerator of agreement
  # 4 - Denominator of agreement
  #Calculates a generalized agreement metric as per France and Carroll (2007)
  #Written Stephen France 2011: france@uwm.edu
  #Please cite
  #France, S.L., Carroll, J.D. (2007).  Development of an Agreement Metric 
  #based upon the RAND Index for the Evaluation of Dimensionality Reduction 
  #Techniques, with Applications to Mapping Customer Data. Machine Learning 
  #and Data Mining in Pattern Recognition, Petra Perner (Ed.), LNAI, 
  #Vol. 4571, Proceedings Conference MLDM 2007, Leipzig/Germany, Springer: Heidelberg
  
  
  n<-nrow(Config1)    #Size for the configuration
  
  #Calculate nearest neighbors
  D<-as.matrix(dist(Config1, method = DMetric, diag = TRUE, upper = TRUE,p=varargin))
  Temp<-apply(D,1,sort,index.return=TRUE)
  IDX1<-t(matrix(unlist(Temp),n*2,n))
  IDX1 <- IDX1[,(n+2):(2*n)]
  
  D<-as.matrix(dist(Config2, method = DMetric, diag = TRUE, upper = TRUE,p=varargin))
  Temp<-apply(D,1,sort,index.return=TRUE)
  IDX2<-t(matrix(unlist(Temp),n*2,n))
  IDX2 <- IDX2[,(n+2):(2*n)]
  rm (D,Temp)
  
  Nok<-Endk-Startk+1
  RBuild<-as.data.frame(matrix(0,Nok,4))
  RowTotalAgree<-matrix(0,n,Nok+1)
  RowAgree<-matrix(0,n,Nok+1)
  TotalAgree<-0
  
  #Update total agreement if not starting at 1
  if (Startk>1)
    {
    for (nCounter in 1:n)
      {
      Temp<-intersect(IDX1[nCounter,1:(Startk-1)],IDX2[nCounter,1:(Startk-1)])
      RowTotalAgree[nCounter,1]<-length(Temp)
      RowAgree[nCounter,1]<-RowTotalAgree[nCounter,1]/(Startk-1)
      TotalAgree<-TotalAgree+RowAgree[nCounter,1]
      }
    rm(Temp)
    }
  
  for (kCounter in Startk:Endk)
  {
    for (nCounter in 1:n)
    {
      #Check for agreement with added item
      Increment<-0
      Increment<-Increment+as.numeric(IDX1[nCounter,kCounter]==IDX2[nCounter,kCounter])
      Increment<-Increment+(length(which(IDX1[nCounter,0:(kCounter-1)]==IDX2[nCounter,kCounter]))>0)
      Increment<-Increment+(length(which(IDX2[nCounter,0:(kCounter-1)]==IDX1[nCounter,kCounter]))>0)
      TotalAgree<-TotalAgree+Increment
      RowTotalAgree[nCounter,kCounter-Startk+2]<-RowTotalAgree[nCounter,kCounter-Startk+1]+Increment
      RowAgree[nCounter,kCounter-Startk+2]<-RowTotalAgree[nCounter,kCounter-Startk+2]/kCounter
    }
    RBuild[kCounter-Startk+1,1]<-TotalAgree/(n*kCounter)
    RBuild[kCounter-Startk+1,2]<-(kCounter)/(n-1)
  }     
      
  #Numerator of aggregate agreement rate
  RBuild[,3]=RBuild[,1]-RBuild[,2]
  #Denominator of aggregate agreement rate 
  RBuild[,4]=1-RBuild[,2]
  RBuild<-cbind(Startk:Endk,RBuild)   
  names(RBuild)<-c('k','Agr','ER(Agr)','AdjNum','AdjDenom')
  
  #Get the final Agreement rate
  OutStruct<-NULL
  OutStruct$n<-n
  OutStruct$Description<-Description
  OutStruct$Startk<-Startk
  OutStruct$Endk<-Endk
  OutStruct$Nok<-Nok
  #Looking for the item names for the items.  
  if (is.null(rownames(Config1))==FALSE) 
  {
    OutStruct$ItemNames<-rownames(Config1)
  } 
  else
  {
    if (is.null(rownames(Config2))==FALSE)
    {
      OutStruct$ItemNames<-rownames(Config2)
    }
  }
  OutStruct$RBuild <-RBuild
  OutStruct$RowTotalAgree<-RowTotalAgree
  OutStruct$RowAgree<-RowAgree
  OutStruct$TotalAgree<-TotalAgree
  OutStruct$Agree<-sum(RBuild[,2])/Nok
  OutStruct$AdjAgree <- sum(RBuild[,4])/sum(RBuild[,5])
  
  return(OutStruct)
}
