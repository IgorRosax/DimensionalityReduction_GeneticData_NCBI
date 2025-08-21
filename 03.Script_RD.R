remove(list = ls())

# Functions

# ************************************************************************************************* #

TimeMeasurementFunction <- function (expr , numberExecutions = 3){
  
  executionTime = rbind()
  i = 0
  while(i < numberExecutions){
    executionTime = rbind(executionTime, system.time(expr()))  
    i = i + 1
  }
  
  return (list('numberOfExecutions' = numberExecutions, 
               'varElapsedTime' = var(executionTime[1:numberExecutions, 3]),
               'meanElapsedTime' = mean(executionTime[1:numberExecutions, 3]),
               'elapsedTimePerExecution' = executionTime[1:numberExecutions, 3]))
}

# ************************************************************************************************* #

gen_DRs <- function  (dirResult, dataSetId){
  load(paste("./datasets/",dataSetId,".RData",sep = ""))
  
  library(dimRed)
  library(smacof)
  library(Rtsne)
  library(umap)
  
  library(smacofx)
  library(RcppHSLocalMDS)
  
  library(pcaMethods)
  
  set.seed(1122)
  
  smacofx_LMDS_parSel <- function(data, conf = NULL, Rn, k, itmax = 10000, 
                                  smallerunitfree = 0.0001, ratio = sqrt(10), n_t = 8)
  {
    
    unitfree<-smallerunitfree*ratio**(0:n_t)
    unitfree<-unitfree[order(unitfree,decreasing = TRUE)]
    
    if(is.null(conf)){
      conf = cmdscale(d = data, k = Rn)
      conf = as.matrix(conf$points)
    }
    
    melhorresul =RcppHSLocalMDS::RcppGetLocalContinuityMetaCriterion(data = data, conf = conf, Rn = Rn, k = k)
    initConf = conf
    BestRD = list("conf" = initConf,
                  "bestTau" = 0,
                  "LocalContinuityResult" = melhorresul)
    for(i in 1:length(unitfree)) {
      resultRD<-smacofx::lmds(delta = data, init = initConf, k = k, tau = unitfree[i], ndim = Rn, itmax = itmax)
      
      CurrentLC<- RcppHSLocalMDS::RcppGetLocalContinuityMetaCriterion(data = data, conf = initConf, Rn = Rn, k = k)
      
      if(CurrentLC$Nk>melhorresul$Nk){
        BestRD = list("conf" = resultRD$conf,
                      "bestTau" = unitfree[i],
                      "LocalContinuityResult" = CurrentLC)
        melhorresul<-CurrentLC
      }}
    return (BestRD)
  }
  
  datasetName = dataSetId
  
  if (!dir.exists(paste("./", dirResult,datasetName, sep = "/"))) dir.create(paste(".",dirResult, datasetName, sep = "/"), recursive = TRUE)
  
  numberExecutions = 1
  
  #dataset preparacao
  
  dataset = list()
  dataset <- as.matrix(dados) #preenche dataset
  dataset <- scale(dataset) # normaliza as dimensões
  datasetDist<- as.matrix(getEuclideanDistanceMatrix(dataset)) #preenche matriz com distancias
  
  
  Rn = 2
  minTau = 0.0001
  ratio = sqrt(10)
  nTau=8
  maxIt = 10000
  optimMethod = "CG"
  num_cpus = 6
  
  
  
  valores_k <- seq(from = Rn + 1, to = (nrow(dataset) - 1), length.out = min((nrow(dataset) - 1) - (Rn ), 30))
  valores_k <- unique(round(valores_k, 0))
  
  conf = cmdscale(datasetDist, Rn)
  conf = as.matrix(conf$conf)
  
  message('Knn Based DRs')
  for (valor_k in valores_k) {
    message('K = ', valor_k)
    
    variacaoKReducaoVsKQuality = 3
    
    for ( kQuality in max(1,valor_k-variacaoKReducaoVsKQuality): min(nrow(datasetDist) - 1, valor_k + variacaoKReducaoVsKQuality))
    {
      # 1. Rcpp Lmds (c++)
      message('1. Rcpp Lmds (c++) [k=',valor_k,']')
      TimeMeasurement = TimeMeasurementFunction( function(){RcppHSlocalMDS(data = datasetDist, conf = as.matrix(conf), Rn = Rn, Kproj = valor_k, Kquality = kQuality, verbose = TRUE,selectBetterUnitFree = TRUE, smallerUnitFree = minTau, n_t = nTau, ratio = ratio, applyHiperbolicSmoothing = FALSE,maxIt = maxIt, optMethod = optimMethod)}, numberExecutions)
      RcppLocalMDSResult = RcppHSlocalMDS(data = datasetDist, conf = as.matrix(conf), Rn = Rn, Kproj = valor_k, Kquality = kQuality, verbose = TRUE,selectBetterUnitFree = TRUE, smallerUnitFree = minTau, n_t = nTau, ratio = ratio, applyHiperbolicSmoothing = FALSE,maxIt = maxIt, optMethod = optimMethod)
      
      #RcppLocalMDSResult$dataset = dataset
      #RcppLocalMDSResult$datasetDist = datasetDist
      RcppLocalMDSResult$TimeMeasurement= TimeMeasurement
      RcppLocalMDSResult$categorias = categorias
      RcppLocalMDSResult$k = valor_k
      RcppLocalMDSResult$kQuality = kQuality
      RcppLocalMDSResult$method = paste("LMDS", " [k=",valor_k,"]", " [kQuality=",kQuality,"]"," (Rcpp)",sep = "")
      
      rdFileName = paste(RcppLocalMDSResult$method,"RData",sep = ".")
      saveRDS(RcppLocalMDSResult, file = paste(".",dirResult,datasetName,rdFileName, sep = "/"))
      rm(rdFileName)
      rm(RcppLocalMDSResult)
      
      # 2. Rcpp HSLmds (c++)
      message('2. Rcpp HSLmds (c++) [k=',valor_k,']')
      
      TimeMeasurement = TimeMeasurementFunction( function(){RcppHSlocalMDS(data = datasetDist, conf = as.matrix(conf), Rn = Rn, Kproj = valor_k, Kquality = kQuality, verbose = TRUE, selectBetterUnitFree = TRUE, smallerUnitFree = minTau, n_t = nTau, ratio = ratio, applyHiperbolicSmoothing = TRUE, gamma = mean(dataset), n_gamma = maxIt, rho = 0.5, maxIt = maxIt, optMethod = optimMethod)}, numberExecutions)
      RcppHSlocalMDSResult = RcppHSlocalMDS(data = datasetDist, conf = as.matrix(conf), Rn = Rn, Kproj = valor_k, Kquality = kQuality, verbose = TRUE, selectBetterUnitFree = TRUE, smallerUnitFree = minTau, n_t = nTau, ratio = ratio, applyHiperbolicSmoothing = TRUE, gamma = mean(dataset), n_gamma = maxIt, rho = 0.5, maxIt = maxIt, optMethod = optimMethod)
      
      #RcppHSlocalMDSResult$dataset = dataset
      #RcppHSlocalMDSResult$datasetDist = datasetDist
      RcppHSlocalMDSResult$TimeMeasurement= TimeMeasurement
      RcppHSlocalMDSResult$categorias = categorias
      RcppHSlocalMDSResult$k = valor_k
      RcppHSlocalMDSResult$kQuality = kQuality
      RcppHSlocalMDSResult$method = paste("HSLMDS", " [k=",valor_k,"]", " [kQuality=",kQuality,"]"," (Rcpp)",sep = "")
      
      rdFileName = paste(RcppHSlocalMDSResult$method,"RData",sep = ".")
      saveRDS(RcppHSlocalMDSResult, file = paste(".",dirResult,datasetName,rdFileName, sep = "/"))
      rm(rdFileName)
      
      rm(RcppHSlocalMDSResult)
      
    }
    
    
    # 3. Smacofx LMDS (R)
    message('3. Smacofx LMDS (R) [k=',valor_k,']')
    TimeMeasurement = TimeMeasurementFunction( function(){smacofx_LMDS_parSel (data = datasetDist, conf = conf, Rn = Rn, k = valor_k, itmax = maxIt, smallerunitfree = minTau, ratio = ratio, n_t = nTau)},numberExecutions)
    smacofxLMDSResult = smacofx_LMDS_parSel (data = datasetDist, conf = conf, Rn = Rn, k = valor_k, itmax = maxIt, smallerunitfree = minTau, ratio = ratio, n_t = nTau)
    
    #smacofxLMDSResult$dataset = dataset
    #smacofxLMDSResult$datasetDist = datasetDist
    smacofxLMDSResultSave <- list()
    smacofxLMDSResultSave$conf = smacofxLMDSResult$conf
    smacofxLMDSResultSave$TimeMeasurement = TimeMeasurement
    smacofxLMDSResultSave$method = paste("LMDS", " [k=",valor_k,"]"," (smacofx R)",sep = "")
    smacofxLMDSResultSave$categorias = categorias
    smacofxLMDSResultSave$k = valor_k
    
    rdFileName = paste(smacofxLMDSResultSave$method,"RData",sep = ".")
    saveRDS(smacofxLMDSResultSave, file = paste(".",dirResult,datasetName,rdFileName, sep = "/"))
    rm(rdFileName)
    rm(smacofxLMDSResult)
    rm(smacofxLMDSResultSave)
    
    
    # 4. UMAP
    message('4. UMAP [k=',valor_k,']')
    TimeMeasurement = TimeMeasurementFunction( function(){umap::umap(d = datasetDist, method = "naive", n_neighbors = valor_k, n_components = Rn)},numberExecutions)
    umapResult <- umap::umap(d = datasetDist, method = "naive", n_neighbors = valor_k, n_components = Rn)
    
    umapResultSet = list()
    #umapResult$dataset = dataset
    #umapResult$datasetDist = datasetDist
    umapResultSet$TimeMeasurement= TimeMeasurement
    umapResultSet$conf = umapResult$layout
    umapResultSet$LocalContinuityResult = RcppGetLocalContinuityMetaCriterion(data = datasetDist, conf = as.matrix(umapResultSet$conf), Rn = Rn, k = valor_k)
    umapResultSet$method = paste("UMAP", " [k=",valor_k,"]",sep = "")
    umapResultSet$categorias = categorias
    umapResultSet$k = valor_k
    
    rdFileName = paste(umapResultSet$method,"RData",sep = ".")
    saveRDS(umapResultSet, file = paste(".",dirResult,datasetName,rdFileName, sep = "/"))
    rm(rdFileName)
    rm(umapResultSet)
    rm(umapResult)
    
    # 5. Isomap
    message('5. Isomap [k=',valor_k,']')
    if( !(datasetName == 'GSE29272') ){
      TimeMeasurement = TimeMeasurementFunction( function(){embed(dataset, "Isomap", .mute = NULL, knn = valor_k, ndim = Rn, get_geod = FALSE)},numberExecutions)
      isomapResult <- embed(dataset, "Isomap", .mute = NULL, knn = valor_k, ndim = Rn, get_geod = FALSE)
      
      isomapResultSet= list()
      #isomapResultSet$object = isomapResult
      #isomapResultSet$dataset = dataset
      isomapResultSet$TimeMeasurement= TimeMeasurement
      isomapResultSet$conf = isomapResult@data@data
      isomapResultSet$LocalContinuityResult = RcppGetLocalContinuityMetaCriterion(data = datasetDist, conf = as.matrix(isomapResultSet$conf), Rn = Rn, k = valor_k)
      isomapResultSet$method = paste("Isomap", " [k=",valor_k,"]",sep = "")
      isomapResultSet$categorias = categorias
      isomapResultSet$k = valor_k
      
      rdFileName = paste(isomapResultSet$method,"RData",sep = ".")
      saveRDS(isomapResultSet, file = paste(".",dirResult,datasetName,rdFileName, sep = "/"))
      rm(rdFileName)
      
      rm(isomapResult)
      rm(isomapResultSet)
    }
    
    # # 6. Hessian Locally Linear Embedding
    # ## Constraints: min(k, n) > d
    # message('6. Hessian Locally Linear Embedding [k=',valor_k,']')
    # 
    # TimeMeasurement = TimeMeasurementFunction( function(){embed(dataset, "HLLE", .mute = NULL, knn = valor_k, ndim = Rn)},numberExecutions)
    # HLLEResult <- embed(dataset, "HLLE", .mute = NULL, knn = valor_k, ndim = Rn)
    # 
    # HLLEResultSet= list()
    # HLLEResultSet$object = HLLEResult
    # HLLEResultSet$dataset = dataset
    # HLLEResultSet$datasetDist = datasetDist
    # HLLEResultSet$TimeMeasurement= TimeMeasurement
    # HLLEResultSet$conf = HLLEResult@data@data
    # HLLEResultSet$LocalContinuityResult = RcppGetLocalContinuityMetaCriterion(data = datasetDist, conf = as.matrix(HLLEResultSet$conf), Rn = Rn, k = valor_k)
    # HLLEResultSet$method = paste("HLLE", " [k=",valor_k,"]",sep = "")
    # HLLEResultSet$categorias = categorias
    # HLLEResultSet$k = valor_k
    # 
    # rdFileName = paste(HLLEResultSet$method,"RData",sep = ".")
    # saveRDS(HLLEResultSet, file = paste(".",dirResult,datasetName,rdFileName, sep = "/"))
    # rm(rdFileName)
    # 
    # rm(HLLEResult)
    # rm(HLLEResultSet)
    
    # 7. t-SNE
    if (valor_k <= 100 & (valor_k<(nrow(dataset)-1)/3)){
      message('7. t-SNE [perp=',valor_k,']')
      
      theta = 0.0
      perp = valor_k
      
      TimeMeasurement = TimeMeasurementFunction( function(){Rtsne::Rtsne(X = dataset, dims = Rn, pca = FALSE, perplexity = perp, theta = theta, num_threads = num_cpus)},numberExecutions)
      tsneResult <- Rtsne::Rtsne(X = dataset, dims = Rn, pca = FALSE, perplexity = perp, theta = theta, num_threads = num_cpus)
      
      tsneResultSet <- list()
      #tsneResultSet$dataset = dataset
      tsneResultSet$TimeMeasurement= TimeMeasurement
      tsneResultSet$conf = tsneResult$Y
      tsneResultSet$LocalContinuityResult = RcppGetLocalContinuityMetaCriterion(data = datasetDist, conf = as.matrix(tsneResultSet$conf), Rn = Rn, k = valor_k)
      tsneResultSet$method = paste("tSNE", " [perp=",valor_k,"]",sep = "")
      tsneResultSet$categorias = categorias
      tsneResultSet$perp = perp
      
      rdFileName = paste(tsneResultSet$method,"RData",sep = ".")
      saveRDS(tsneResultSet, file = paste(".",dirResult,datasetName,rdFileName, sep = "/"))
      rm(rdFileName)
      rm(tsneResult)
      rm(tsneResultSet)
      rm(theta)
      rm(perp)
    }
  }
  
  message('Global Based DRs')
  
  numberExecutions = 3
  k = 20
  
  # 8. PCA SVD
  message('8. PCA SVD')
  
  TimeMeasurement = TimeMeasurementFunction( function(){pcaMethods::pca(dataset, method = "svd", nPcs = Rn)},numberExecutions)
  PcaSvdResult <- pcaMethods::pca(dataset, method = "svd", nPcs = Rn)
  
  PcaSvdResultSet= list()
  #PcaSvdResultSet$dataset = dataset
  PcaSvdResultSet$TimeMeasurement= TimeMeasurement
  PcaSvdResultSet$conf = PcaSvdResult@scores
  PcaSvdResultSet$LocalContinuityResult = RcppGetLocalContinuityMetaCriterion(data = datasetDist, conf = as.matrix(PcaSvdResultSet$conf), Rn = Rn, k = k)
  PcaSvdResultSet$method = "PCA_SVD"
  PcaSvdResultSet$categorias = categorias
  
  saveRDS(PcaSvdResultSet, file = paste(".",dirResult,datasetName,"PCA_SVD.RData", sep = "/"))
  rm(PcaSvdResult)
  rm(PcaSvdResultSet)
  
  # 9. PPCA
  message('9. PPCA')
  
  TimeMeasurement = TimeMeasurementFunction( function(){pcaMethods::pca(dataset, method = "ppca", nPcs = Rn)},numberExecutions)
  PPcaResult <- pcaMethods::pca(dataset, method = "ppca", nPcs = Rn)
  
  PPcaResultSet= list()
  #PPcaResultSet$dataset = dataset
  PPcaResultSet$TimeMeasurement= TimeMeasurement
  PPcaResultSet$conf = PPcaResult@scores
  PPcaResultSet$LocalContinuityResult = RcppGetLocalContinuityMetaCriterion(data = datasetDist, conf = as.matrix(PPcaResultSet$conf), Rn = Rn, k = k)
  PPcaResultSet$method = "PPCA"
  PPcaResultSet$categorias = categorias
  
  saveRDS(PPcaResultSet, file = paste(".",dirResult,datasetName,"PPCA.RData", sep = "/"))
  rm(PPcaResult)
  rm(PPcaResultSet)
  
  
  # 10. Kernel PCA
  message('10. Kernel PCA')
  TimeMeasurement = TimeMeasurementFunction( function(){embed(dataset, "kPCA", .mute = NULL, knn = 5, ndim = 1)},numberExecutions)
  kPcaResult <- embed(dataset, "kPCA", .mute = NULL, knn = 5, ndim = 1)
  
  kPcaResultSet= list()
  #kPcaResultSet$dataset = dataset
  kPcaResultSet$TimeMeasurement= TimeMeasurement
  kPcaResultSet$conf = kPcaResult@data@data
  kPcaResultSet$LocalContinuityResult = RcppGetLocalContinuityMetaCriterion(data = datasetDist, conf = as.matrix(kPcaResultSet$conf), Rn = Rn, k = k)
  kPcaResultSet$method = "kPCA"
  kPcaResultSet$categorias = categorias
  
  saveRDS(kPcaResultSet, file = paste(".",dirResult,datasetName,"kPCA.RData", sep = "/"))
  rm(kPcaResultSet)
  rm(kPcaResult)
  
  
  # 11. Diffusion Maps 
  message('11. Diffusion Maps ')
  TimeMeasurement = TimeMeasurementFunction( function(){embed(dataset, "DiffusionMaps", .mute = NULL, ndim = Rn)},numberExecutions)
  DiffusionMapsResult <- embed(dataset, "DiffusionMaps", .mute = NULL, ndim = Rn)
  
  DiffusionMapsResultSet= list()
  #DiffusionMapsResultSet$dataset = dataset
  DiffusionMapsResultSet$TimeMeasurement= TimeMeasurement
  DiffusionMapsResultSet$conf = DiffusionMapsResult@data@data
  DiffusionMapsResultSet$LocalContinuityResult = RcppGetLocalContinuityMetaCriterion(data = datasetDist, conf = as.matrix(DiffusionMapsResultSet$conf), Rn = Rn, k = k)
  DiffusionMapsResultSet$method = "DiffusionMaps"
  DiffusionMapsResultSet$categorias = categorias
  
  saveRDS(DiffusionMapsResultSet, file = paste(".",dirResult,datasetName,"DiffusionMaps.RData", sep = "/"))
  rm(DiffusionMapsResult)
  rm(DiffusionMapsResultSet)
  
  # 12. Dimensionality Reduction via Regression  
  message('12. Dimensionality Reduction via Regression  ')
  
  TimeMeasurement = TimeMeasurementFunction( function(){embed(dataset, "DRR", .mute = NULL, ndim = Rn)}, numberExecutions)
  DdrResult <- embed(dataset, "DRR", .mute = NULL, ndim = Rn)
  
  DdrResultSet= list()
  #DdrResultSet$dataset = dataset
  DdrResultSet$TimeMeasurement= TimeMeasurement
  DdrResultSet$conf = DdrResult@data@data
  DdrResultSet$LocalContinuityResult = RcppGetLocalContinuityMetaCriterion(data = datasetDist, conf = as.matrix(DdrResultSet$conf), Rn = Rn, k = k)
  DdrResultSet$method = "DRR"
  DdrResultSet$categorias = categorias
  
  saveRDS(DdrResultSet, file = paste(".",dirResult,datasetName,"DRR.RData", sep = "/"))
  rm(DdrResult)
  rm(DdrResultSet)
  
  
  # 13. Multidimensional Scaling Smacof
  message('13. Multidimensional Scaling Smacof')
  TimeMeasurement = TimeMeasurementFunction( function(){smacof::mds(delta = datasetDist, ndim = Rn, type = "ratio")},numberExecutions)
  smacofMDSResult <- smacof::mds(delta = datasetDist, ndim = Rn, type = "ratio")
  
  #smacofMDSResult$dataset = dataset
  smacofMDSResult$TimeMeasurement= TimeMeasurement
  smacofMDSResult$LocalContinuityResult = RcppGetLocalContinuityMetaCriterion(data = datasetDist, conf = as.matrix(smacofMDSResult$conf), Rn = Rn, k = k)
  smacofMDSResult$method = "smacof_MDS (R)"
  smacofMDSResult$categorias = categorias
  
  saveRDS(smacofMDSResult, file = paste(".",dirResult,datasetName,"MDS smacof (R).RData", sep = "/"))
  rm(smacofMDSResult)
  
  # 14. MDS
  message('14. MDS')
  TimeMeasurement = TimeMeasurementFunction( function(){RcppHSMDS(data = datasetDist, conf = conf, Rn = Rn, Kquality = k, verbose = TRUE, applyHiperbolicSmoothing = FALSE, maxIt = maxIt, optMethod = "CG")},numberExecutions)
  RcppKruskalMDSResult = RcppHSMDS(data = datasetDist, conf = conf, Rn = Rn, Kquality = k, verbose = TRUE, applyHiperbolicSmoothing = FALSE, maxIt = maxIt, optMethod = "CG")
  
  #RcppKruskalMDSResult$dataset = dataset
  #RcppKruskalMDSResult$datasetDist = datasetDist
  RcppKruskalMDSResult$TimeMeasurement = TimeMeasurement
  RcppKruskalMDSResult$method = "MDS (Rcpp)"
  RcppKruskalMDSResult$categorias = categorias
  
  saveRDS(RcppKruskalMDSResult, file = paste(".",dirResult,datasetName,"MDS (Rcpp).RData", sep = "/"))
  rm(RcppKruskalMDSResult)
  
  # 15. HS MDS
  message('15. HS MDS')
  TimeMeasurement = TimeMeasurementFunction( function(){RcppHSMDS(data = datasetDist, conf = conf, Rn = Rn, Kquality = k, verbose = TRUE, applyHiperbolicSmoothing = TRUE, gamma = mean(datasetDist), n_gamma = maxIt, rho = 0.5, maxIt = maxIt, optMethod = "CG")},numberExecutions)
  RcppHSMDSResult = RcppHSMDS(data = datasetDist, conf = conf, Rn = Rn, Kquality = k, verbose = TRUE, applyHiperbolicSmoothing = TRUE, gamma = mean(datasetDist), n_gamma = maxIt, rho = 0.5, maxIt = maxIt, optMethod = "CG")
  
  #RcppHSMDSResult$dataset = dataset
  #RcppHSMDSResult$datasetDist = datasetDist
  RcppHSMDSResult$TimeMeasurement = TimeMeasurement
  RcppHSMDSResult$method = "HSMDS (Rcpp)"
  RcppHSMDSResult$categorias = categorias
  
  saveRDS(RcppHSMDSResult, file = paste(".",dirResult,datasetName,"HSMDS (Rcpp).RData", sep = "/"))
  rm(RcppHSMDSResult)
}

# ************************************************************************************************* #

# Main

series = rbind('GSE14020','GSE28735','GSE18842','GSE35988','GSE21034','GSE44076','GSE29272')
for (serie in series){
  
  dataId = serie
  
  unlink(paste(".","RD_Results",dataId,"*",sep = "/"), recursive = TRUE)
  
  message('DataSet: ', serie)
  
  gen_DRs(dirResult = "RD_Results", dataSetId =  dataId)
  
}

# ************************************************************************************************* #
