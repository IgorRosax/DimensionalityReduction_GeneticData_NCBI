PLiftSinglegg <- function(AgreeStruct,ShowMode="show",IsFill=TRUE,TextSize=10){
  #ScatterSinglegg
  #Creates a lift diagram which shows the agreement of a technique over k 
  #(the k defined by the agreement structures) in comparison to expected agreement
  #AgreeStruct - An agreement structure created with GenAgree or GenAgreeDist
  #IsFill - TRUE if the areas between the agreement curves and 
  # the expected agreement are colored (with transparent colors for overlapping areas)
  #         FALSE otherwise
  #ShowMode - "show" if wish to show plot in graphics environment
  #          - "return" if wish to return the plot 
  #TextSize - The text size in standard ggplot units
  #Written Stephen France 2016-2019: sfrance@business.msstat.edu
  #Please cite
  #France, S.L., Carroll, J.D. (2007).  Development of an Agreement Metric 
  #based upon the RAND Index for the Evaluation of Dimensionality Reduction 
  #Techniques, with Applications to Mapping Customer Data. Machine Learning 
  #and Data Mining in Pattern Recognition, Petra Perner (Ed.), LNAI, 
  #Vol. 4571, Proceedings Conference MLDM 2007, Leipzig/Germany, Springer: Heidelberg  
  
  require(ggplot2)
  require(reshape2)
  
  #Build a data frame, with agreement for the configuration and expected agreement
  PlotVal<-AgreeStruct$RBuild[1:3]  
  PlotTitle<-paste("Agree:",round(AgreeStruct$Agree,4)," AdjAgree:",round(AgreeStruct$AdjAgree,4))
  PlotVal <- cbind(PlotVal,min_line=pmin(PlotVal[,2],PlotVal[,3]) ) 
  PlotVal <- melt(PlotVal, id.vars=c("k","min_line"), variable.name="Config", value.name="Agree")
  PlotVal$Config<-as.factor(PlotVal$Config)
  
  if (IsFill==TRUE)
  {
    #Create filled plot
    sp <- ggplot(data=PlotVal, aes(x=k, fill=Config,color=Config,linetype=Config))
    sp <- sp + geom_ribbon(aes(ymax=Agree, ymin=min_line))
    sp <- sp + geom_line(aes(y = Agree),size=1)
    sp <- sp + scale_linetype_manual(values=c(Agr="solid", `ER(Agr)`="twodash"))
    sp <- sp + scale_fill_manual(values=c(Agr="darkblue", `ER(Agr)`="white"))
    sp <- sp + scale_color_manual(values=c(Agr="black", `ER(Agr)`="darkred"))
  }
  else
  {
    #Create a standard line plot
    sp <- ggplot(data=PlotVal, aes(x=k,color=Config,linetype=Config))
    sp <- sp + geom_line(aes(y = Agree),size=1)
    sp <- sp + scale_color_manual(values=c(Agr="blue",`ER(Agr)`="darkred"))
    sp <- sp + scale_linetype_manual(values=c(Agr="solid", `ER(Agr)`="twodash"))
  }
  
  sp <- sp + ggtitle(PlotTitle)
  sp <-sp + theme_bw()
  sp <-sp + theme(legend.title=element_text(face="bold"))
  sp <-sp + theme(text=element_text(size = TextSize),plot.title = element_text(hjust = 0.5,size = rel(1),face="bold"))

  if (ShowMode=="show")
  {
    plot(sp)
  }
  else
  {
    return(sp)
  }  
  
  
}