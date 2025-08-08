

remove(list = ls())

library("GEOquery")
seriesMat = rbind(cbind(id = 'GSE14020', clsVar = "source_name_ch1", hasNaAttr = FALSE),
                  cbind(id = 'GSE28735', clsVar = "characteristics_ch1", hasNaAttr = FALSE),
                  cbind(id = 'GSE18842', clsVar = "characteristics_ch1.1", hasNaAttr = FALSE),
                  cbind(id = 'GSE35988', clsVar = "characteristics_ch2", hasNaAttr = TRUE),
                  cbind(id = 'GSE21034', clsVar = "characteristics_ch1.3", hasNaAttr = FALSE),
                  cbind(id = 'GSE44076', clsVar = "source_name_ch1", hasNaAttr = FALSE),
                  cbind(id = 'GSE29272', clsVar = "source_name_ch1", hasNaAttr = FALSE),
                  cbind(id = 'GSE39582', clsVar = "source_name_ch1", hasNaAttr = FALSE))

for (row in (1:nrow(seriesMat))){
  serie <- seriesMat[row,]
  load(paste("datasets/",serie["id"],".RData",sep = ""))
  # variavel_classes <- serie["clsVar"]
  # 
  print( paste("Dataset: ", serie["id"], "(",nrow(as.data.frame(dados)),"x",ncol(as.data.frame(dados)),")", sep = "" ) )
  if (serie["hasNaAttr"])
    print("TRUE")
  print( "Classes: "  )
  classes <- summary(factor(categorias))
  print (classes)
  #summary(factor(pData(phenoData(matriz_informacoes[[1]]))[[variavel_classes]]))
  print("")
  print("")
  # summary(factor(pData(phenoData(matriz_informacoes[[1]]))[[variavel_classes]]))
  
}
