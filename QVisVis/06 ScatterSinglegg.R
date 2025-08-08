ScatterSinglegg <- function(AgreeStruct,Dest,PlotMode="simple",LowCol="darkred",MidCol="white",HighCol="darkblue",
                            Startk=NULL,Endk=NULL,ShowMode="show",GradientMode="none",TextSize=10,TitleElements=c(TRUE,FALSE,FALSE,TRUE,TRUE)){
  #ScatterSinglegg
  #Plots the agreement results at a disaggregate level.  Plots results in a 2D scatterplot
  #for each value of i and k
  #AgreeStruct - An agreement structure created with GenAgree or GenAgreeDist
  #Dest - The destination mapping
  #PlotMode - The mode for the heatmap
  # "simple" - A simple plot, showing agreement from 0 to 1 in shades of a single color
  # "compare" - Subtract the expected value of agreement for each value of k
  # "adjusted" - The adjusted agreement, (Agr-E(Agr))/(1-E(Agr)). Is infinity for k-1
  # LowCol,MidCol,HighCol - The plot colors. Default to red(-ve), white(0), and blue(+ve)
  # Startk,Endk - Optional k values to allow the plot to be restricted
  # ShowMode - "show" if wish to show plot in graphics environment
  #          - "return" if wish to return the plot 
  # GradientMode - "none" - Just show the points
  #              - "loess" - Show a loess transformation on the points
  # TextSize - The text size in standard ggplot units
  # TitleElements - Array of 5 items. 1st technique name, 2nd dataset name, 3rd k-values, 4th agreement, 5th adjusted agreement
  #Written Stephen France 2016: sfrance@business.msstat.edu
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
  
  PlotTitle<-paste(kText,"Agree:",round(PlAgree,4)," AdjAgree:",round(PlAdjAgree,4))
  
  Dest<-as.data.frame(Dest)
  names(Dest)<-c("D1","D2")
  
  PlotVal<-AgreeStruct$RowAgree[,(Startk+1):(Endk+1),drop=FALSE]
  if ((PlotMode=="compare")||(PlotMode=="adjusted"))
  {
    ExpVal<-outer(rep(1,n),((Startk:Endk)/(n-1)))
    if (PlotMode=="compare")
    {
      PlotVal<-(PlotVal-ExpVal)
      AgreeName<-"IncAgree"
    }
    else
    {
      PlotVal<-(PlotVal-ExpVal)/(1-ExpVal)
      AgreeName<-"AdjAgree"
      #This will give infinity values, which will mess up the plot
      if (Endk==(n-1))
      {
        PlotVal=PlotVal[,1:AgreeStruct$Nok-1]
      }
      
    }
  }
  else
  {
    AgreeName<-"Agreement"
  }
  #In this case, we need an aggregate row per item
  PlotVal<-apply(PlotVal,1,mean,drop=FALSE)
  Dest<-cbind(Dest,PlotVal)
  names(Dest)[3]<-"Agree"
  
  if (GradientMode=="none")
  {
    sp<-ggplot(data =  Dest, aes(x = D1, y = D2,color=Agree))+geom_point()
    sp<-sp+scale_color_gradient2(low = LowCol,mid=MidCol, high=HighCol,guide="colorbar",name=AgreeName,midpoint=0)
    sp<-sp+ggtitle(PlotTitle)
    sp<-sp+theme(legend.title=element_text(face="bold"))
    sp<-sp+theme_bw()
    sp<-sp+theme(legend.title=element_text(face="bold"))
    sp<-sp+theme(text=element_text(size = TextSize),plot.title = element_text(hjust = 0.5,size = rel(1),face="bold"))
  }
  else
  {
    #create the loess scaling output
    loess.out <- loess(Agree ~ D1 * D2, data = Dest)
    RangeD1<-(max(Dest$D1)-min(Dest$D1))
    RangeD2<-(max(Dest$D2)-min(Dest$D2))
    xgrid <-  seq(min(Dest$D1), max(Dest$D1),RangeD1/50)
    ygrid <-  seq(min(Dest$D2), max(Dest$D2),RangeD2/50)
    
    NewData <-  expand.grid(D1 = xgrid, D2 = ygrid)
    ShadeAgree<-  predict(loess.out, newdata = NewData)
    
    #Create an interpolation for the shaded data
    ShadeAgree2 <- melt(ShadeAgree, id.vars = c("D1", "D2"), measure.vars = Agree)
    ShadeAgree2$D1<-as.numeric(substr(ShadeAgree2$D1, 4,20))
    ShadeAgree2$D2<-as.numeric(substr(ShadeAgree2$D2, 4,20))
    colnames(ShadeAgree2)[3]<-"Agree"
    
    #Create the ggplot output
    sp <- ggplot(ShadeAgree2, aes(x = D1, y = D2, z = Agree))
    sp <-sp+stat_contour(geom = "polygon", aes(fill = ..level..))
    sp <-sp+geom_tile(aes(fill = Agree)) 
    sp <-sp+stat_contour(bins = 15)
    sp <-sp+scale_fill_gradient2(low = LowCol,mid=MidCol, high=HighCol,guide="colorbar",name="Agreement\n",midpoint=0)
    sp <-sp+guides(fill = guide_colorbar(title = "Agree"))
    sp <-sp+geom_point(data = Dest, aes(x = D1, y = D2,fill=Agree), shape = 21, size = 2,color="black")
    sp <-sp+scale_color_gradient2(low = LowCol,mid=MidCol, high=HighCol,midpoint=0)
    sp<-sp+ggtitle(PlotTitle)
    sp<-sp+theme_bw()
    sp<-sp+theme(legend.title=element_text(face="bold"))
    sp<-sp+theme(text=element_text(size = TextSize),plot.title = element_text(hjust = 0.5,size = rel(1),face="bold"))
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