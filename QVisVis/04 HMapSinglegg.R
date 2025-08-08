HMapSinglegg <- function(AgreeStruct,PlotMode="simple",LowCol="darkred",MidCol="white",HighCol="darkblue",
                         Startk=NULL,Endk=NULL,ShowMode="show",NoxLabels=10,NoyLabels=10,
                         TitleElements=c(TRUE,FALSE,FALSE,TRUE,TRUE),Ordering=NULL,TextSize=10){
  #HMapSinglegg
  #Plots the agreement results at a disaggregate level.  Plots results in a heatmap
  #for each value of i and k
  #AgreeStruct - An agreement structure created with GenAgree or GenAgreeDist
  #PlotMode - The mode for the heatmap
  # "simple" - A simple plot, showing agreement from 0 to 1 in shades of a single color
  # "compare" - Subtract the expected value of agreement for each value of k
  # "adjusted" - The adjusted agreement, (Agr-E(Agr))/(1-E(Agr)). Is infinity for k-1
  # LowCol,MidCol,HighCol - The plot colors. Default to red(-ve), white(0), and blue(+ve)
  # Startk,Endk - Optional k values to allow the plot to be restricted
  # ShowMode - TRUE if show graph and FALSE if return graph as function output
  # TitleElements - Array of 5 items. 1st technique name, 2nd dataset name, 3rd k-values, 4th agreement, 5th adjusted agreement
  # Ordering - A double or integer ordering, the same length as the number of items. Defaults to NULL.
  # TextSize - The text size in standard ggplot units
  #Written Stephen France 2016-2019: sfrance@business.msstat.edu
  #Please cite
  #France, S.L., Carroll, J.D. (2007).  Development of an Agreement Metric 
  #based upon the RAND Index for the Evaluation of Dimensionality Reduction 
  #Techniques, with Applications to Mapping Customer Data. Machine Learning 
  #and Data Mining in Pattern Recognition, Petra Perner (Ed.), LNAI, 
  #Vol. 4571, Proceedings Conference MLDM 2007, Leipzig/Germany, Springer: Heidelberg
  require(ggplot2)
  require(reshape2)
  
  if (is.null(Startk)&&is.null(Endk))
    {
    Fixk=FALSE
    kText=''
    }
  else {
    Fixk=TRUE}
  if (TitleElements[1]==TRUE)
  {
    kText=AgreeStruct$Description
  }
  else
  {
    kText=""
  }
  if (TitleElements[2]==TRUE)
  {
    kText=AgreeStruct$DataDescription
  }
  
  n<-AgreeStruct$n
  if (is.null(Startk))
    {Startk<-AgreeStruct$Startk}
  if (is.null(Endk)) 
    {Endk<-AgreeStruct$Endk}
  if (Fixk==TRUE)
  { 
    if (TitleElements[3]==TRUE)
    {
    if (Startk==Endk)
      {kAdd=paste("k=",Startk," ",sep="")}
    else
      {kAdd=paste("k=",Startk,":",Endk," ",sep="")}
      kText=paste(k,kAdd,sep=" ")
    }
    AgreeStruct<-SubAgreeDist(AgreeStruct,Startk,Endk)
    PlAgree<-AgreeStruct$SubAgree
    PlAdjAgree<-AgreeStruct$SubAdjAgree
  }
  else
  {
    PlAgree<-AgreeStruct$Agree
    PlAdjAgree<-AgreeStruct$AdjAgree
  }
  if (TitleElements[4]==TRUE)
  {
    PlotTitle<-paste(kText,"Agree:",round(PlAgree,4),sep=" ")
  }
  if (TitleElements[5]==TRUE)
  {
    PlotTitle<-paste(kText,"AdjAgree:",round(PlAdjAgree,5),sep=" ")
  }
  
  #Row agree goes base row then Startk to Endk
  PlotVal<-AgreeStruct$RowAgree[,(Startk+1):(Endk+1),drop=FALSE]
  if ((PlotMode=="compare")||(PlotMode=="adjusted"))
  {
    ExpVal<-outer(rep(1,n),((Startk:Endk)/(n-1)))
    if (PlotMode=="compare")
    {
      PlotVal<-(PlotVal-ExpVal)
    }
    else
    {
      PlotVal<-(PlotVal-ExpVal)/(1-ExpVal)
    }
  }
  
  #If ordered then reorder plotval
  if (is.null(Ordering)==FALSE)
  {
    PlotVal<-PlotVal[Ordering,]
  }
  
  PlotVal<-as.data.frame(cbind(1:n,PlotVal))
  
  names(PlotVal)<-c("Item",as.character(Startk:Endk))
  PlotVal$Item<-as.factor(PlotVal$Item)
  
  #Create a structure with colummn for item,  column for k, and a column for the row agreement
  PlotVal <- melt(PlotVal, id.vars=c("Item"), variable.name="k", value.name="RowAgree")
  PlotVal$k<-as.factor(PlotVal$k)
  
  sp<-ggplot(data =  PlotVal, aes(x = k, y = Item))
  sp<-sp+geom_tile(aes(fill = RowAgree))
  sp<-sp+theme(legend.title=element_text(face="bold"))
  sp<-sp+scale_fill_gradient2(low = LowCol,mid=MidCol, high=HighCol, guide="colorbar",name="Agreement\n",midpoint=0)
  sp<-sp+ggtitle(PlotTitle)
  sp<-sp+theme_bw()
  sp<-sp+theme(legend.title=element_text(face="bold"))
  sp<-sp+theme(text=element_text(size = TextSize),plot.title = element_text(hjust = 0.5,size = rel(1),face="bold"))
  
  if (n>NoyLabels)
  {
    Jump<-round(n/NoyLabels)
    sp<-sp+scale_y_discrete(breaks=seq(0,n,by=Jump))
  }
  if (AgreeStruct$Nok>NoxLabels)
  {
    Jump<-round(AgreeStruct$Nok/NoxLabels)
    sp<-sp+scale_x_discrete(breaks=seq(Startk,Endk,by=Jump))
  }
  if (ShowMode=="show")
  {
    plot(sp)
  }
  else
  {
    return(sp)
  }
}