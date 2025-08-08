# GEOquery.R
# Este arquivo baixa os bancos de dados disponíveis no GEO.

#Este arquivo foi modificado para baixar os conjuntos de dados desejados listados em 'series'
remove(list = ls())
library("GEOquery")

series = rbind(cbind(id = 'GSE14020', clsVar = "source_name_ch1", hasNaAttr = FALSE),
               cbind(id = 'GSE28735', clsVar = "characteristics_ch1", hasNaAttr = FALSE),
               cbind(id = 'GSE18842', clsVar = "characteristics_ch1.1", hasNaAttr = FALSE),
               cbind(id = 'GSE35988', clsVar = "characteristics_ch2", hasNaAttr = TRUE),
               cbind(id = 'GSE21034', clsVar = "characteristics_ch1.3", hasNaAttr = FALSE),
               cbind(id = 'GSE44076', clsVar = "source_name_ch1", hasNaAttr = FALSE),
               cbind(id = 'GSE29272', clsVar = "source_name_ch1", hasNaAttr = FALSE),
               cbind(id = 'GSE39582', clsVar = "source_name_ch1", hasNaAttr = FALSE))

for (row in (1:nrow(series))){
  serie <- series[row,]
  
  diretorio_BDs <- "datasets"
  if (!dir.exists(diretorio_BDs)) dir.create(diretorio_BDs, recursive = TRUE)
  arquivo_saida <- paste(serie["id"], "RData", sep = ".")
  path_saida <- paste(".", diretorio_BDs, arquivo_saida, sep = "/")
  
  matriz_informacoes <- getGEO(serie["id"], GSEMatrix = TRUE)
  
  nomes_amostras <- pData(phenoData(matriz_informacoes[[1]]))$geo_accession
  gsm_amostras <- lapply(nomes_amostras, getGEO, GSEMatrix = TRUE)
  
  numero_de_genes <- as.numeric(phenoData(matriz_informacoes[[1]])$data_row_count[1])
  
  dados <- rep(0, numero_de_genes)
  for (i in seq_along(nomes_amostras)) {
    dados <- cbind(dados, Table(gsm_amostras[i][[1]])[2])
  }
  dados <- t(dados[-1])
  
  if (serie["hasNaAttr"]){
    dados <- as.data.frame(dados)
    dados <- dados[, colSums(is.na(dados)) == 0]
  }
  
  categorias <- pData(phenoData(matriz_informacoes[[1]]))[[serie["clsVar"]]]
  
  save.image(file = path_saida)
}


