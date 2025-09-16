remove(list = ls())

library(HSLocalMDS)
library(stringr)
library(dplyr)
library(RColorBrewer) 
diretorio_fonte_QVisVis <- "./QVisVis"
GenAgree_source <- paste(diretorio_fonte_QVisVis, "03 GenAgree.R", sep = "/")
HMapSinglegg_source <- paste(diretorio_fonte_QVisVis, "04 HMapSinglegg.R", sep = "/")

source(GenAgree_source, encoding = "UTF8")
source(HMapSinglegg_source, encoding = "UTF8")

color_palette <- c(
  brewer.pal(8, "Dark2"), 
  brewer.pal(7, "Set1")
)

shape_palette <- rep(c(15, 16, 17, 18, 19), times = 3)
# ************************************************************************************************* #
# Functions
# ************************************************************************************************* #

run_Knn <- function ( Conf, ResultsetClasses, MethodName ){
  
  library(caret)
  
  if (nlevels(factor(ResultsetClasses)) < 2){
    print(paste("knn not eligible. there must be at least 2 factors levels in the data. [",MethodName,"].",sep = ""))
    return()
  }
  
  set.seed(1234)
  
  Categ = as.data.frame(ResultsetClasses)
  colnames(Categ) = "Categorias"
  Categ$Categorias[Categ$Categorias==""] <- "NA"
  
  classifData = cbind(Conf,Categ)
  
  
  intrain <- createDataPartition(y = factor(classifData[,"Categorias"]), p= 0.7, list = FALSE)
  training <- classifData[intrain,]
  testing <- classifData[-intrain,]
  
  #dim(training); dim(testing);
  
  knn_fit <- caret::train( Categorias ~ ., data = training, method = "knn",
                           trControl=trainControl(method = "repeatedcv", number = 10, repeats = 3, classProbs = FALSE, allowParallel = TRUE),
                           tuneGrid = expand.grid(k = seq(1,max(summary(factor(training$Categorias))),by = 3)))
  
  knn_predict <- predict(knn_fit, testing[,1:ncol(testing)-1])
  
  confusionMat <- confusionMatrix(knn_predict, factor(testing$Categorias))
  
  return (data.frame(  Accuracy    = round(confusionMat$overall["Accuracy"],3),
                       Sensitivity = round(confusionMat$byClass["Sensitivity"],3),
                       Specificity = round(confusionMat$byClass["Specificity"],3),
                       Precision   = round(confusionMat$byClass["Precision"],3),
                       Recall      = round(confusionMat$byClass["Recall"],3),
                       F1          = round(confusionMat$byClass["F1"],3),
                       stringsAsFactors = FALSE,
                       row.names = MethodName) )
  
  
}

# ************************************************************************************************* #

bestLCMC <- function (dataDist, conf){
  
  Rn = 2
  
  k_start <- 1
  k_end <- nrow(dataDist) - 1
  
  for ( k in k_start:k_end){
    LCMC =HSLocalMDS::RcppGetLocalContinuityMetaCriterion(data = dataDist, conf = conf, Rn = Rn, k = k)
    if (k == 1L){
      BestLCMC <- LCMC
      BestLCMC$k = k
    }
    if(LCMC$Mk_adjusted > BestLCMC$Mk_adjusted  ){
      BestLCMC <- LCMC
      BestLCMC$k = k
    }
  }
  
  return (data.frame(  best_K_evaluation    = BestLCMC$k,
                       LCMC_Nk = round(BestLCMC$Nk,3),
                       LCMC_Mk = round(BestLCMC$Mk,3),
                       LCMC_Mk_Adjusted   = round(BestLCMC$Mk_adjusted,3),
                       stringsAsFactors = FALSE) )
}

# ************************************************************************************************* #

genAgreeAvgComparisonGraphic <- function (resultsTable,diretorio_output, datasetName){
  
  knnBasedResults <- resultsTable %>%
    mutate(k = coalesce(k_reduction, perp))%>%
    filter(!is.na(k)) %>%
    group_by(method) %>%
    filter(n_distinct(k) > 1) %>%
    ungroup() %>%
    mutate(tipo_linha = "Knn Based")
  
  min_k <- min(knnBasedResults$k, na.rm = TRUE)
  max_k <- max(knnBasedResults$k, na.rm = TRUE)
  
  notKnnBasedResults <- resultsTable %>%
    mutate(k = coalesce(k_reduction, perp))%>%
    group_by(method) %>%
    filter(all(is.na(k))) %>%
    reframe(
      AgreeAvg = first(AgreeAvg), 
      k = c(min_k, max_k)
    ) %>%
    ungroup() %>%
    mutate(tipo_linha = "Not Knn Based")
  
  resultsBinded <- bind_rows(
    select(knnBasedResults, method, k, AgreeAvg, tipo_linha),
    select(notKnnBasedResults, method, k, AgreeAvg, tipo_linha)
  )
  
  final_plot <- ggplot(
    data = resultsBinded, 
    aes(x = k, y = AgreeAvg, color = method, shape = method, group = method, linetype = tipo_linha)
  ) +
    geom_line(size = 1) +
    geom_point(data = knnBasedResults, size = 2.5) +
    scale_linetype_manual(values = c("Knn Based" = "solid", "Not Knn Based" = "dashed")) +
    labs(
      #title = "Comparação de Métodos de Redução de Dimensionalidade",
      
      x = "Neighborhood/Perplexity (k)",
      y = "Average Agreement Rate (AgreeAvg)",
      color = "Method",
      shape = "Method",
      linetype = "Method type"
    ) +
    theme_minimal() +
    scale_color_manual(values = color_palette) +
    scale_shape_manual(values = shape_palette) +
    theme(
      legend.position = "bottom",
      legend.box = "vertical",
      #plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
      #plot.subtitle = element_text(hjust = 0.5, size = 12)
    )
  fileNamePlot = paste("07.agreeAvgComparison_k_perp_",datasetName,".png",sep="")
  ggsave(paste(diretorio_output, datasetName,fileNamePlot,sep = "/"), plot = final_plot, width = 7, height = 4)
}

# ************************************************************************************************* #

evaluate_DRs <- function (dirResult, datasetName){
  load(paste("./datasets/",datasetName,".RData",sep = ""))
  dataset <- as.matrix(dados) #preenche dataset
  dataset <- scale(dataset) # normaliza as dimensões
  datasetDist<- as.matrix(getEuclideanDistanceMatrix(dataset)) #preenche matriz com distancias

  diretorio_output = "./output"
  
  if (!dir.exists(paste(diretorio_output, datasetName, sep = "/"))) dir.create(paste(diretorio_output, datasetName, sep = "/"), recursive = TRUE)
  
  arquivos_RDs <- list.files(path = paste(".", dirResult,datasetName, sep = "/"), pattern = ".*\\.RData$", full.names = TRUE)  
  
  resultsTable <- data.frame(matrix(ncol = 19, nrow = 0))  
  colnames(resultsTable) <- c("id", "method","k_reduction", "perp", "AgreeAvg","Phi_Fx","numberOfExecutions", "varElapsedTime", "meanElapsedTime", "Accuracy", "Sensitivity", "Specificity", "Precision", "Recall", "F1","best_K_evaluation","LCMC_Nk","LCMC_Mk","LCMC_Mk_Adjusted")
  
  for (arquivo in arquivos_RDs) {
    message('Evaluating DR: ', arquivo)
    
    path_arquivo <- paste(".", arquivo, sep = "/")
    obj <- readRDS(path_arquivo)
    if (ncol(obj$conf) < 2){
      obj$conf <- cbind(obj$conf, rep(0, length(obj$conf)))
    }
    
    if (!is.null(obj$k)){
      k_reduction = obj$k
    }else{
      k_reduction = NA
    }
    
    if (!is.null(obj$perp)){
      perp = obj$perp
    }else{
      perp = NA
    }
    
    methodNameClean = gsub(" \\[.*\\]", "", obj$method)
    
    k_lim_inferior <- 1
    k_lim_superior <- nrow(dataset) - 1
    
    genAgreeResult <- GenAgree(Config1 = dataset, Config2 = obj$conf, Startk = k_lim_inferior, Endk = k_lim_superior, DMetric = "euclidian")
    DR_results <- data.frame(  id = obj$method,
                               method = methodNameClean,
                               k_reduction = k_reduction,
                               perp = perp,
                               AgreeAvg = round(genAgreeResult$Agree,2),
                               Phi_Fx = round(genAgreeResult$AdjAgree,2),
                               numberOfExecutions = obj$TimeMeasurement$numberOfExecutions,
                               varElapsedTime = round(obj$TimeMeasurement$varElapsedTime,3),
                               meanElapsedTime = round(obj$TimeMeasurement$meanElapsedTime,3),
                               stringsAsFactors = FALSE,
                               row.names = obj$method)
    
    knnResult = run_Knn( Conf = obj$conf, ResultsetClasses = obj$categorias, obj$method)
    
    if(is.null(datasetDist))
      bestLcmcResult = bestLCMC(dataDist = as.matrix(getEuclideanDistanceMatrix(dataset)), conf = obj$conf)
    else
      bestLcmcResult = bestLCMC(dataDist = datasetDist, conf = obj$conf)
    
    
    resultsTable <- rbind( resultsTable,
                           cbind(DR_results,
                                 knnResult,
                                 bestLcmcResult) )
    
    ## Gerando Graficos
    
    # # gerando o scatterplot com os dados em dimensão reduzida e com cores de acordo com classificação dos dados.
    # nome_arquivo_scatterplot_classified <- paste(obj$method, "scatterplot_classified")
    # nome_arquivo_scatterplot_classified <- paste(nome_arquivo_scatterplot_classified, "png", sep = ".")
    # nome_arquivo_scatterplot_classified <- paste(diretorio_output, datasetName, nome_arquivo_scatterplot_classified, sep = "/")
    # png(filename = nome_arquivo_scatterplot_classified)
    # plot(obj$conf[,1], obj$conf[,2], col = factor(obj$categorias) ,asp=1,cex.sub=0.7,cex.main=0.7, xlab="x", ylab="y")
    # dev.off()
    
    # # gerando o heatmap
    # heatmap <- HMapSinglegg(AgreeStruct = genAgreeResult, ShowMode = "return")
    # nome_arquivo_heatmap <- paste(obj$method, "heatmap")
    # nome_arquivo_heatmap <- paste(nome_arquivo_heatmap, "png", sep = ".")
    # nome_arquivo_heatmap <- paste(diretorio_output, datasetName, nome_arquivo_heatmap, sep = "/")
    # ggsave(nome_arquivo_heatmap, plot = heatmap)
    
    # gerando curva agree
    agreeCurveResult <- cbind(id= obj$method ,genAgreeResult$RBuild)
    
    if (exists("agreeCurveResultsBinded")){
      agreeCurveResultsBinded <- rbind(agreeCurveResultsBinded,agreeCurveResult)
    }else{
      agreeCurveResultsBinded <- agreeCurveResult
    }
    
    # agreeCurve <- ggplot(
    #   data = agreeCurveResult,
    #   aes(x = k, y = Agr)
    # ) +
    #   geom_line(size = 1) +
    #   geom_point(data = agreeCurveResult, size = 2.5) +
    #   labs(
    #     title = "Agreement Rate Curve",
    #     x = "Neighborhood (k)",
    #     y = "Agreement Rate"
    #   ) +
    #   theme_minimal() +
    #   scale_color_manual(values = color_palette) +
    #   scale_shape_manual(values = shape_palette) +
    #   theme(
    #     legend.position = "bottom",
    #     legend.box = "vertical",
    #     plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    #     plot.subtitle = element_text(hjust = 0.5, size = 12)
    #   )
    # fileNamePlot = paste(obj$method,"agreeCurve.png",sep = "")
    # ggsave(paste(diretorio_output, datasetName,fileNamePlot,sep = "/"), plot = agreeCurve, width = 11, height = 8)
    
  }
  
  arquivo_tabela <- paste("resultsTable_", datasetName, ".csv", sep = "")
  arquivo_tabela <- paste(diretorio_output, datasetName, arquivo_tabela, sep = "/")
  write.csv(resultsTable, file = arquivo_tabela)
  
  
  
  # #best_AgreeAvg
  # best_records <- resultsTable %>%
  #   group_by(method) %>%
  #   arrange(desc(AgreeAvg), desc(Phi_Fx), desc(Accuracy), meanElapsedTime) %>%
  #   slice_head(n = 1) %>%
  #   ungroup() %>%
  #   arrange(desc(AgreeAvg), desc(Phi_Fx), desc(Accuracy), meanElapsedTime) 
  # 
  # arquivo_tabela <- paste("resultsTable_best_AgreeAvg_", datasetName, ".csv", sep = "")
  # arquivo_tabela <- paste(diretorio_output, datasetName, arquivo_tabela, sep = "/")
  # write.csv(best_records, file = arquivo_tabela)
  # 
  # comparativeAgreeCurve <- ggplot(
  #   data = semi_join(agreeCurveResultsBinded, best_records, by = "id"),
  #   aes(x = k, y = Agr, color = id, shape = id, group = id)
  # ) +
  #   geom_line(size = 1) +
  #   geom_point(data = semi_join(agreeCurveResultsBinded, best_records, by = "id"), size = 2.5) +
  #   labs(
  #     title = "Agreement Rate Curve",
  #     x = "Neighborhood (k)",
  #     y = "Agreement Rate",
  #     color = "Method",
  #     shape = "Method"
  #   ) +
  #   theme_minimal() +
  #   scale_color_manual(values = color_palette) +
  #   scale_shape_manual(values = shape_palette) +
  #   theme(
  #     legend.position = "bottom",
  #     legend.box = "vertical",
  #     plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
  #     plot.subtitle = element_text(hjust = 0.5, size = 12)
  #   )
  # fileNamePlot = paste("01.comparativeAgreeCurve_best_AgreeAvg_",datasetName,".png",sep="")
  # ggsave(paste(diretorio_output, datasetName,fileNamePlot,sep = "/"), plot = comparativeAgreeCurve, width = 11, height = 8)
  
  # #best_Phi_Fx
  # best_records <- resultsTable %>%
  #   group_by(method) %>%
  #   arrange(desc(Phi_Fx), desc(AgreeAvg), desc(Accuracy), meanElapsedTime) %>%
  #   slice_head(n = 1) %>%
  #   ungroup() %>%
  #   arrange(desc(Phi_Fx), desc(AgreeAvg), desc(Accuracy), meanElapsedTime) 
  # 
  # arquivo_tabela <- paste("resultsTable_best_Phi_Fx_", datasetName, ".csv", sep = "")
  # arquivo_tabela <- paste(diretorio_output, datasetName, arquivo_tabela, sep = "/")
  # write.csv(best_records, file = arquivo_tabela)
  # 
  # comparativeAgreeCurve <- ggplot(
  #   data = semi_join(agreeCurveResultsBinded, best_records, by = "id"),
  #   aes(x = k, y = Agr, color = id, shape = id, group = id)
  # ) +
  #   geom_line(size = 1) +
  #   geom_point(data = semi_join(agreeCurveResultsBinded, best_records, by = "id"), size = 2.5) +
  #   labs(
  #     title = "Agreement Rate Curve",
  #     x = "Neighborhood (k)",
  #     y = "Agreement Rate",
  #     color = "Method",
  #     shape = "Method"
  #   ) +
  #   theme_minimal() +
  #   scale_color_manual(values = color_palette) +
  #   scale_shape_manual(values = shape_palette) +
  #   theme(
  #     legend.position = "bottom",
  #     legend.box = "vertical",
  #     plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
  #     plot.subtitle = element_text(hjust = 0.5, size = 12)
  #   )
  # fileNamePlot = paste("02.comparativeAgreeCurve_best_Phi_Fx_",datasetName,".png",sep="")
  # ggsave(paste(diretorio_output, datasetName,fileNamePlot,sep = "/"), plot = comparativeAgreeCurve, width = 11, height = 8)
  
  # #best_Accuracy
  # best_records <- resultsTable %>%
  #   group_by(method) %>%
  #   arrange(desc(Accuracy), desc(AgreeAvg), desc(Phi_Fx), meanElapsedTime) %>%
  #   slice_head(n = 1) %>%
  #   ungroup() %>%
  #   arrange(desc(Accuracy),desc(AgreeAvg), desc(Phi_Fx), meanElapsedTime) 
  # 
  # arquivo_tabela <- paste("resultsTable_best_Accuracy_", datasetName, ".csv", sep = "")
  # arquivo_tabela <- paste(diretorio_output, datasetName, arquivo_tabela, sep = "/")
  # write.csv(best_records, file = arquivo_tabela)
  # 
  # comparativeAgreeCurve <- ggplot(
  #   data = semi_join(agreeCurveResultsBinded, best_records, by = "id"),
  #   aes(x = k, y = Agr, color = id, shape = id, group = id)
  # ) +
  #   geom_line(size = 1) +
  #   geom_point(data = semi_join(agreeCurveResultsBinded, best_records, by = "id"), size = 2.5) +
  #   labs(
  #     title = "Agreement Rate Curve",
  #     x = "Neighborhood (k)",
  #     y = "Agreement Rate",
  #     color = "Method",
  #     shape = "Method"
  #   ) +
  #   theme_minimal() +
  #   scale_color_manual(values = color_palette) +
  #   scale_shape_manual(values = shape_palette) +
  #   theme(
  #     legend.position = "bottom",
  #     legend.box = "vertical",
  #     plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
  #     plot.subtitle = element_text(hjust = 0.5, size = 12)
  #   )
  # fileNamePlot = paste("03.comparativeAgreeCurve_best_Accuracy_",datasetName,".png",sep="")
  # ggsave(paste(diretorio_output, datasetName,fileNamePlot,sep = "/"), plot = comparativeAgreeCurve, width = 11, height = 8)
  
  # #best_LCMC
  # best_records <- resultsTable %>%
  #   group_by(method) %>%
  #   arrange(desc(LCMC_Mk_Adjusted),desc(LCMC_Mk), desc(AgreeAvg), desc(Accuracy), meanElapsedTime) %>%
  #   slice_head(n = 1) %>%
  #   ungroup() %>%
  #   arrange(desc(LCMC_Mk_Adjusted),desc(LCMC_Mk), desc(AgreeAvg), desc(Accuracy), meanElapsedTime) 
  # 
  # arquivo_tabela <- paste("resultsTable_best_LCMC_", datasetName, ".csv", sep = "")
  # arquivo_tabela <- paste(diretorio_output, datasetName, arquivo_tabela, sep = "/")
  # write.csv(best_records, file = arquivo_tabela)
  # 
  # comparativeAgreeCurve <- ggplot(
  #   data = semi_join(agreeCurveResultsBinded, best_records, by = "id"),
  #   aes(x = k, y = Agr, color = id, shape = id, group = id)
  # ) +
  #   geom_line(size = 1) +
  #   geom_point(data = semi_join(agreeCurveResultsBinded, best_records, by = "id"), size = 2.5) +
  #   labs(
  #     title = "Agreement Rate Curve",
  #     x = "Neighborhood (k)",
  #     y = "Agreement Rate",
  #     color = "Method",
  #     shape = "Method"
  #   ) +
  #   theme_minimal() +
  #   scale_color_manual(values = color_palette) +
  #   scale_shape_manual(values = shape_palette) +
  #   theme(
  #     legend.position = "bottom",
  #     legend.box = "vertical",
  #     plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
  #     plot.subtitle = element_text(hjust = 0.5, size = 12)
  #   )
  # fileNamePlot = paste("04.comparativeAgreeCurve_best_LCMC_",datasetName,".png",sep="")
  # ggsave(paste(diretorio_output, datasetName,fileNamePlot,sep = "/"), plot = comparativeAgreeCurve, width = 11, height = 8)
  
  # #best_agreeRate_variousKRed
  # best_records <- inner_join(agreeCurveResultsBinded, resultsTable, by = "id") %>%
  #   group_by(method,k) %>%
  #   arrange(desc(Agr), desc(`ER(Agr)`)) %>%
  #   slice_head(n = 1) %>%
  #   ungroup() %>%
  #   arrange(method,k)
  # 
  # comparativeAgreeCurve <- ggplot(
  #   data = best_records,
  #   aes(x = k, y = Agr, color = method, shape = method, group = method)
  # ) +
  #   geom_line(size = 1) +
  #   geom_point(data = best_records, size = 2.5) +
  #   labs(
  #     title = paste("Agreement Rate Curve (", datasetName,")",sep=""),
  #     x = "Neighborhood (k)",
  #     y = "Agreement Rate (AR)",
  #     color = "Method",
  #     shape = "Method"
  #   ) +
  #   theme_minimal() +
  #   scale_color_manual(values = color_palette) +
  #   scale_shape_manual(values = shape_palette) +
  #   theme(
  #     legend.position = "bottom",
  #     legend.box = "vertical",
  #     plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
  #     plot.subtitle = element_text(hjust = 0.5, size = 12)
  #   )
  # fileNamePlot = paste("05.comparativeAgreeCurve_best_Agree_kRed_",datasetName,".png",sep="")
  # ggsave(paste(diretorio_output, datasetName,fileNamePlot,sep = "/"), plot = comparativeAgreeCurve, width = 11, height = 8)
  
  #best_agreeRateAdjusted_variousKRed
  best_records <- inner_join(agreeCurveResultsBinded, resultsTable, by = "id") %>%
    mutate(AgrAdj = Agr - `ER(Agr)`) %>%
    group_by(method,k) %>%
    arrange(desc(AgrAdj), desc(Agr), desc(`ER(Agr)`)) %>%
    slice_head(n = 1) %>%
    ungroup() %>%
    arrange(method,k)
  
  comparativeAgreeCurve <- ggplot(
    data = best_records,
    aes(x = k, y = AgrAdj, color = method, shape = method, group = method)
  ) +
    geom_line(size = 1) +
    geom_point(data = best_records, size = 2.5) +
    labs(
      #title = paste("Agreement Rate Adjusted Curve (", datasetName,")",sep=""),
      x = "Neighborhood (k)",
      y = "Agreement Rate Adjusted (AR*)",
      color = "Method",
      shape = "Method"
    ) +
    theme_minimal() +
    scale_color_manual(values = color_palette) +
    scale_shape_manual(values = shape_palette) +
    theme(
      legend.position = "bottom",
      legend.box = "vertical",
      #plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
      #plot.subtitle = element_text(hjust = 0.5, size = 12)
    )
  fileNamePlot = paste("06.comparativeAgreeCurve_best_AgreeAdjusted_kRed_",datasetName,".png",sep="")
  ggsave(paste(diretorio_output, datasetName,fileNamePlot,sep = "/"), plot = comparativeAgreeCurve, width = 7, height = 4)
  
  #comparisonGraphics
  genAgreeAvgComparisonGraphic(resultsTable = resultsTable, diretorio_output = diretorio_output, datasetName = datasetName)
}


# ************************************************************************************************* #
# Main
# ************************************************************************************************* #
dirResult = "RD_Results"


series = rbind('GSE14020','GSE28735','GSE18842','GSE35988','GSE21034','GSE44076','GSE29272')
for (serie in series){
  
  dataId = serie
  
  unlink(paste(".","output",dataId,"*",sep = "/"), recursive = TRUE)
  
  message('DataSet: ', serie)
  
  #evaluate_DRs
  evaluate_DRs(dirResult, dataId)
  
}