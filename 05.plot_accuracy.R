remove(list = ls())

library(NLDR)
library(stringr)
library(dplyr)
library(RColorBrewer) 
library(ggplot2)
library(tidyr) # Para a função pivot_longer, equivalente ao melt do Python
library(reshape2)

diretorio_fonte_QVisVis <- "./QVisVis"
GenAgree_source <- paste(diretorio_fonte_QVisVis, "03 GenAgree.R", sep = "/")
HMapSinglegg_source <- paste(diretorio_fonte_QVisVis, "04 HMapSinglegg.R", sep = "/")

source(GenAgree_source, encoding = "UTF8")
source(HMapSinglegg_source, encoding = "UTF8")

color_palette <- c(
  brewer.pal(8, "Dark2"), 
  brewer.pal(8, "Set1")
)

shape_palette <- rep(c(0,1,2,5,6,15,16,17,18), times = 2)

#methodsList <- levels(factor(resultsTable$method))
methodsList <-c("DiffusionMaps", "DRR", "HSLMDS (NLDR)", "HSMDS (NLDR)", "Isomap", "kPCA", "LMDS (NLDR)", "LMDS (smacofx)", 
                "MDS (NLDR)", "MDS (smacof)", "PCA_SVD", "PPCA","tSNE", "UMAP")#, "HSLMDS (rNLDR)", "LMDS (rNLDR)")

color_palette <- setNames(color_palette, methodsList)
shape_palette <- setNames(shape_palette, methodsList)

methodsListRenamed <-c("DiffusionMaps" = "DiffusionMaps", 
                       "DRR" = "DRR", 
                       "HSLMDS (NLDR)" = "HSLocalMDS (NLDR)", 
                       #"HSLMDS (rNLDR)" = "HSLocalMDS (R)", 
                       "HSMDS (NLDR)" = "HSMDS (NLDR)", 
                       "Isomap" = "Isomap", 
                       "kPCA" = "kPCA",
                       "LMDS (NLDR)" = "LocalMDS (NLDR)",
                       #"LMDS (rNLDR)" = "LocalMDS (R)",
                       "LMDS (smacofx)" = "LocalMDS (smacofx)", 
                       "MDS (NLDR)" = "MDS (NLDR)", 
                       "MDS (smacof)" = "MDS (smacof)", 
                       "PCA_SVD" = "PCA", 
                       "PPCA" = "PPCA",
                       "tSNE" = "t-SNE", 
                       "UMAP" = "UMAP")

diretorio_output = "./output"

resultsTable <- data.frame()
series = rbind('GSE14020','GSE28735','GSE18842','GSE35988','GSE21034','GSE44076','GSE29272')
for (datasetName in series){
  
  datasetName
  
  message('DataSet: ', datasetName)
  
  fileNameCSV = paste("06.comparativeAgreeAdjusted_kRed_",datasetName,".csv",sep="")
  AgreeAdjResults = read.csv(paste(diretorio_output,datasetName,fileNameCSV, sep = "/"), sep = "," )
  
  comparativeAgreeCurve <- ggplot(
    data = subset(AgreeAdjResults, method %in% names(methodsListRenamed)),
    aes(x = k, y = AgrAdj, color = method, shape = method, group = method)
  ) +
    geom_line(linewidth = 0.5) +
    geom_point(data = . %>% filter(k %% max(1, round(max(AgreeAdjResults$k) / 50)) == 0 | k == min(k) | k == max(k)), 
               size = 2) +
    labs(
      #title = paste("Agreement Rate Adjusted Curve (", datasetName,")",sep=""),
      x = "Vizinhança (k)",
      y = "Taxa de Concordância Ajustada (AR*)",
      color = "Método",
      shape = "Método"
    ) +
    theme_minimal() +
    scale_color_manual(
      values = color_palette,
      labels = methodsListRenamed) +
    scale_shape_manual(
      values = shape_palette,
      labels = methodsListRenamed) +
    guides(
      color = guide_legend(nrow = 16, ncol = 1),
      shape = guide_legend(nrow = 16, ncol = 1)) +
    theme(
      legend.position = "right",
      legend.box = "vertical",
      #plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
      #plot.subtitle = element_text(hjust = 0.5, size = 12)
    )
  fileNamePlot = paste("06.comparativeAgreeCurve_best_AgreeAdjusted_kRed_",datasetName,".png",sep="")
  ggsave(paste(diretorio_output, datasetName,fileNamePlot,sep = "/"), plot = comparativeAgreeCurve, width = 10, height = 6)
  
  rm(fileNameCSV)
  rm(AgreeAdjResults)
  rm(comparativeAgreeCurve)
  rm(fileNamePlot)
  
  load(paste("./datasets/",datasetName,".RData",sep = ""))
  
  resultsTable =rbind( resultsTable , 
                       cbind( dataset = datasetName,
                              n_samples = nrow(dados),
                              n_attr = ncol(dados),
                              n_values = nrow(dados) * ncol(dados),
                              read.csv(paste(diretorio_output,"/",datasetName,"/resultsTable_",datasetName,".csv", sep = ""), sep = "," ) ) )
  
  
}


timeResults <- resultsTable %>%
  group_by(dataset,method,n_samples, n_attr, n_values) %>%
  summarise(
    meanElapsedTime = mean(meanElapsedTime, trim = 0.01, na.rm = TRUE),
    .groups = 'drop' 
  ) %>%
  arrange(dataset)


performance_plot <- ggplot(
  data = subset(timeResults, method %in% names(methodsListRenamed)), 
  aes(x = n_samples, y = meanElapsedTime, color = method, shape = method, group = method)
) +
  geom_line(linewidth = 0.5) +
  geom_point(size = 2) +
  
  #facet_wrap(~ dataset, scales = "free_x") + 
  labs(
    #title = "Dimentionality Reduction Methods Performance by Sample Size",
    #subtitle = "Mean execution time comparison",
    x = "Número de Amostras (n)",
    y = "Média do Tempo de Execução (seconds)",
    color = "Método",
    shape = "Método"
  ) +
  theme_minimal() +
  scale_color_manual(
    values = color_palette,
    labels = methodsListRenamed) +
  scale_shape_manual(
    values = shape_palette,
    labels = methodsListRenamed) +
  guides(
    color = guide_legend(nrow = 16, ncol = 1),
    shape = guide_legend(nrow = 16, ncol = 1)) +
  theme(
    legend.position = "right",
    legend.box = "vertical",
    #plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    #plot.subtitle = element_text(hjust = 0.5, size = 12)
  )
# Exiba o gráfico

ggsave(paste(diretorio_output,"10.timeComparison_method_n.png",sep = "/"), plot = performance_plot, width = 10, height = 6)


performance_plot <- ggplot(
  data =  subset(timeResults, method %in% names(methodsListRenamed)), 
  aes(x = n_values, y = meanElapsedTime, color = method, shape = method, group = method)
) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 2.5) +
  
  #facet_wrap(~ dataset, scales = "free_x") + 
  labs(
    title = "Dimentionality Reduction Methods Performance by Dataset Size",
    subtitle = "Mean execution time comparison",
    x = "Number of Samples * Number of genes",
    y = "Mean Execution Time (seconds)",
    color = "Method",
    shape = "Method"
  ) +
  theme_minimal() +
  scale_color_manual(
    values = color_palette,
    labels = methodsListRenamed) +
  scale_shape_manual(
    values = shape_palette,
    labels = methodsListRenamed) +
  theme(
    legend.position = "bottom",
    legend.box = "vertical",
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    plot.subtitle = element_text(hjust = 0.5, size = 12)
  )

# Exiba o gráfico

ggsave(paste(diretorio_output,"11.timeComparison_method_nXgen.png",sep = "/"), plot = performance_plot, width = 7, height = 4)




#Knn Results

#best_Accuracy
best_records <- resultsTable %>%
  group_by(dataset,method) %>%
  arrange(desc(Accuracy), desc(AgreeAvg), desc(Phi_Fx), meanElapsedTime) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  arrange(dataset, desc(Accuracy),desc(AgreeAvg), desc(Phi_Fx), meanElapsedTime)

write.csv(best_records, file = paste(diretorio_output, "01.resultsTable_best_Accuracy.csv", sep = "/"))


# #plot_lines
# best_Accuracy_Plot <- ggplot(
#   data = best_records,
#   aes(x = dataset, y = Accuracy, color = method, shape = method, linetype = method, group = method)
# ) +
#   geom_line(size = 1) +
#   geom_point(data = best_records, size = 2.5) +
#   labs(
#     title = "Agreement Rate Curve",
#     x = "Neighborhood (k)",
#     y = "Agreement Rate",
#     color = "Method",
#     shape = "Method",
#     linetype = "Method"
#   ) +
#   theme_minimal() +
#   scale_color_manual(
#     values = color_palette,
#     labels = methodsListRenamed) +
#   scale_shape_manual(
#     values = shape_palette,
#     labels = methodsListRenamed) +
#   theme(
#     legend.position = "bottom",
#     legend.box = "vertical",
#     plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
#     plot.subtitle = element_text(hjust = 0.5, size = 12)
#   )
# fileNamePlot = paste("03.comparativeAgreeCurve_best_Accuracy_",datasetName,".png",sep="")
# ggsave(paste(diretorio_output, datasetName,fileNamePlot,sep = "/"), plot = comparativeAgreeCurve, width = 7, height = 4)
# 
# 
# plot(best_Accuracy_Plot)

# Plot piles

#Accuracy
df_melted_acc <- best_records %>%
  select(dataset, method, Accuracy) %>%
  pivot_longer(
    cols = Accuracy,
    names_to = "Metrica",
    values_to = "Valor"
  )

chart_acc <- ggplot(
  subset(df_melted_acc, method %in% names(methodsListRenamed)),
  aes(x = method, y = Valor, fill = dataset)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(
    #title = "Média(AR_k) por Método e Dataset",
    x = "Método",
    y = "Acurácia",
    fill = "Dataset"
  ) +
  theme_minimal() +
  scale_color_manual(values = color_palette) +
  scale_shape_manual(values = shape_palette) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + # Rotaciona rótulos do eixo X
  scale_x_discrete(labels = methodsListRenamed) +
  scale_y_continuous(labels = scales::comma) # Formata o eixo Y para números grandes

ggsave(paste(diretorio_output, "01.comparative_best_Accuracy.png",sep = "/"), plot = chart_acc, width = 5.5, height = 3.5)

#Sensitivity

best_records <- resultsTable %>%
  group_by(dataset,method) %>%
  arrange(desc(Sensitivity), desc(AgreeAvg), desc(Phi_Fx), meanElapsedTime) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  arrange(dataset, desc(Sensitivity),desc(AgreeAvg), desc(Phi_Fx), meanElapsedTime)

write.csv(best_records, file = paste(diretorio_output, "02.resultsTable_best_Sensitivity.csv", sep = "/"))

df_melted_Sensitivity <- best_records %>%
  select(dataset, method, Sensitivity) %>%
  pivot_longer(
    cols = Sensitivity,
    names_to = "Metrica",
    values_to = "Valor"
  )

chart_Sensitivity <- ggplot(
  subset(df_melted_Sensitivity, method %in% names(methodsListRenamed)),
  aes(x = method, y = Valor, fill = dataset)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(
    #title = "Média(AR_k) por Método e Dataset",
    x = "Método",
    y = "Sensibilidade",
    fill = "Dataset"
  ) +
  theme_minimal() +
  scale_color_manual(values = color_palette) +
  scale_shape_manual(values = shape_palette) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + # Rotaciona rótulos do eixo X
  scale_x_discrete(labels = methodsListRenamed) +
  scale_y_continuous(labels = scales::comma) # Formata o eixo Y para números grandes

ggsave(paste(diretorio_output, "02.comparative_best_Sensitivity.png",sep = "/"), plot = chart_Sensitivity, width = 5.5, height = 3.5)

#Specificity
best_records <- resultsTable %>%
  group_by(dataset,method) %>%
  arrange(desc(Specificity), desc(AgreeAvg), desc(Phi_Fx), meanElapsedTime) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  arrange(dataset, desc(Specificity),desc(AgreeAvg), desc(Phi_Fx), meanElapsedTime)

write.csv(best_records, file = paste(diretorio_output, "03.resultsTable_best_Specificity.csv", sep = "/"))

df_melted_Specificity <- best_records %>%
  select(dataset, method, Specificity) %>%
  pivot_longer(
    cols = Specificity,
    names_to = "Metrica",
    values_to = "Valor"
  )

chart_Specificity <- ggplot(
  subset(df_melted_Specificity, method %in% names(methodsListRenamed)),
  aes(x = method, y = Valor, fill = dataset)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(
    #title = "Média(AR_k) por Método e Dataset",
    x = "Método",
    y = "Especificidade",
    fill = "Dataset"
  ) +
  theme_minimal() +
  scale_color_manual(values = color_palette) +
  scale_shape_manual(values = shape_palette) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + # Rotaciona rótulos do eixo X
  scale_x_discrete(labels = methodsListRenamed) +
  scale_y_continuous(labels = scales::comma) # Formata o eixo Y para números grandes

ggsave(paste(diretorio_output, "03.comparative_best_Specificity.png",sep = "/"), plot = chart_Specificity, width = 5.5, height = 3.5)

#Precision
best_records <- resultsTable %>%
  group_by(dataset,method) %>%
  arrange(desc(Precision), desc(AgreeAvg), desc(Phi_Fx), meanElapsedTime) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  arrange(dataset, desc(Precision),desc(AgreeAvg), desc(Phi_Fx), meanElapsedTime)

write.csv(best_records, file = paste(diretorio_output, "04.resultsTable_best_Precision.csv", sep = "/"))

df_melted_Precision <- best_records %>%
  select(dataset, method, Precision) %>%
  pivot_longer(
    cols = Precision,
    names_to = "Metrica",
    values_to = "Valor"
  )

chart_Precision <- ggplot(
  subset(df_melted_Precision, method %in% names(methodsListRenamed)),
  aes(x = method, y = Valor, fill = dataset)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(
    #title = "Média(AR_k) por Método e Dataset",
    x = "Método",
    y = "Precisão",
    fill = "Dataset"
  ) +
  theme_minimal() +
  scale_color_manual(values = color_palette) +
  scale_shape_manual(values = shape_palette) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + # Rotaciona rótulos do eixo X
  scale_x_discrete(labels = methodsListRenamed) +
  scale_y_continuous(labels = scales::comma) # Formata o eixo Y para números grandes

ggsave(paste(diretorio_output, "04.comparative_best_Precision.png",sep = "/"), plot = chart_Precision, width = 5.5, height = 3.5)

#Recall
best_records <- resultsTable %>%
  group_by(dataset,method) %>%
  arrange(desc(Recall), desc(AgreeAvg), desc(Phi_Fx), meanElapsedTime) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  arrange(dataset, desc(Recall),desc(AgreeAvg), desc(Phi_Fx), meanElapsedTime)

write.csv(best_records, file = paste(diretorio_output, "05.resultsTable_best_Recall.csv", sep = "/"))

df_melted_Recall <- best_records %>%
  select(dataset, method, Recall) %>%
  pivot_longer(
    cols = Recall,
    names_to = "Metrica",
    values_to = "Valor"
  )

chart_Recall <- ggplot(
  subset(df_melted_Recall, method %in% names(methodsListRenamed)),
  aes(x = method, y = Valor, fill = dataset)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(
    #title = "Média(AR_k) por Método e Dataset",
    x = "Método",
    y = "Recall",
    fill = "Dataset"
  ) +
  theme_minimal() +
  scale_color_manual(values = color_palette) +
  scale_shape_manual(values = shape_palette) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + # Rotaciona rótulos do eixo X
  scale_x_discrete(labels = methodsListRenamed) +
  scale_y_continuous(labels = scales::comma) # Formata o eixo Y para números grandes

ggsave(paste(diretorio_output, "05.comparative_best_Recall.png",sep = "/"), plot = chart_Recall, width = 5.5, height = 3.5)

#F1
best_records <- resultsTable %>%
  group_by(dataset,method) %>%
  arrange(desc(F1), desc(AgreeAvg), desc(Phi_Fx), meanElapsedTime) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  arrange(dataset, desc(F1),desc(AgreeAvg), desc(Phi_Fx), meanElapsedTime)

write.csv(best_records, file = paste(diretorio_output, "06.resultsTable_best_F1.csv", sep = "/"))

df_melted_F1 <- best_records %>%
  select(dataset, method, F1) %>%
  pivot_longer(
    cols = F1,
    names_to = "Metrica",
    values_to = "Valor"
  )

chart_F1 <- ggplot(
  subset(df_melted_F1, method %in% names(methodsListRenamed)),
  aes(x = method, y = Valor, fill = dataset)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(
    #title = "Média(AR_k) por Método e Dataset",
    x = "Método",
    y = "F1-Score",
    fill = "Dataset"
  ) +
  theme_minimal() +
  scale_color_manual(values = color_palette) +
  scale_shape_manual(values = shape_palette) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + # Rotaciona rótulos do eixo X
  scale_x_discrete(labels = methodsListRenamed) +
  scale_y_continuous(labels = scales::comma) # Formata o eixo Y para números grandes

ggsave(paste(diretorio_output, "06.comparative_best_F1.png",sep = "/"), plot = chart_F1, width = 5.5, height = 3.5)

