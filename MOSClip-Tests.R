## 1. loading the required libraries

library(org.Hs.eg.db)
library(pheatmap)
library(RColorBrewer)
library(EDASeq)
library(graphite)
library(houseOfClipUtility)
library(MOSClip)
library(curatedTCGAData)
library(TCGAutils)
library(tidyverse)
library(dplyr)
library(Hmisc)
library(readxl)
library(data.table)
library(gridExtra)
library(ggVennDiagram)
library(ggplot2)
library(reshape2)
library(igraph)
library(openxlsx)
library(penalized)
library(gridExtra)
library(cowplot)
library(survival)
library(survminer)

## 2. MOSClip pathway test

### 2.1 BRCA

#Filtering the Reactome pathways
#choosing only the ones who  have  expression data available 
nodesLength <- sapply(reactome, function(g) {
  length(intersect(graphite::nodes(g), row.names(BRCA_EXP)))
  })
reactHuge_BRCA <- reactome[nodesLength >= 10]

#Creating the multiomics
multiOmics_BRCA_pathway <- Omics(data = list(expr = BRCA_EXP, met = as.matrix(BRCA_Methylation), mut = BRCA_Mutation, cnv = BRCA_CNV), methods = c("summarizeWithPca", "summarizeInCluster", "summarizeToBinaryEvents", "summarizeToBinaryEvents"), specificArgs = list(pcaArgs = list(name = "exp", shrink = "TRUE", method = "topological", maxPCs = 3), clusterArgs = list(name = "met"), binaryArgs = list(name = "mut", binaryClassMin = 10), binaryArgs = list(name = "cnv", binaryClassMin = 10)))

#Considering only the genes which are present in expression data frame
genesToConsider_BRCA <- row.names(BRCA_EXP)

#Survival pathway test
startTime <- Sys.time()
multiOmicsFullHuge_BRCA <- lapply(reactHuge_BRCA, function(g) {
        print(g@title)
        set.seed(1234)
        fcl = multiOmicsSurvivalPathwayTest(multiOmics_BRCA_pathway, g, survival_BRCA, useThisGenes =genesToConsider_BRCA)
        fcl
    })
endTime <- Sys.time()
print(endTime - startTime)
save(multiOmicsFullHuge_BRCA, file = "./multiOmicsFullHuge_BRCA.RData")

#Testing the result
multiOmicsFullHuge_BRCA <- multiOmicsFullHuge_BRCA[!duplicated(multiOmicsFullHuge_BRCA)]
multiPathwayReport(multiOmicsFullHuge_BRCA)

#Resampling step
pathwaySummaryHuge <- multiPathwayReport(multiOmicsFullHuge_BRCA)
useThisPathwaysFM <- row.names(pathwaySummaryHuge[pathwaySummaryHuge$pvalue <= 0.05, ])
sPathwayM <- pathwaySummaryHuge[useThisPathwaysFM, , drop = T]
startTime <- Sys.time()
permsPathwaysHuge_BRCA <- MOSClip:::resamplingPathway(fullMultiOmics = multiOmics_BRCA_pathway, survival_BRCA, reactHuge_BRCA, nperm = 10, pathwaySubset = useThisPathwaysFM)
endTime <- Sys.time()
print(endTime - startTime)
save(permsPathwaysHuge_BRCA, file = "./permsPathwaysHuge_BRCA.RData")

#Adding the resampling score to the previous matrix
pathwaySummaryHuge <- multiPathwayReport(multiOmicsFullHuge_BRCA)
useThisPathwaysFM <- row.names(pathwaySummaryHuge[pathwaySummaryHuge$pvalue <= 0.05, ])
sPathwayM <- pathwaySummaryHuge[useThisPathwaysFM, , drop = T]
resampligSuccessCount <- getPathwaysModulesSuccess(perms = permsPathwaysHuge_BRCA, moduleSummary = sPathwayM)
pathwaySummaryHuge_BRCA <- addResamplingCounts(pathwaySummaryHuge, resampligSuccessCount)

#Chossing those pathways with a score of at least 8
pathwaySummaryHuge_BRCA <- pathwaySummaryHuge_BRCA[pathwaySummaryHuge_BRCA$resamplingCount >= 8,]
save(pathwaySummaryHuge_BRCA, file = "./pathwaySummaryHuge_BRCA.RData")
pathwaySummaryHuge_BRCA

### 2.2 CESC 

#Filtering the reactome pathways
#choosing only the ones who  have  expression data available
nodesLength <- sapply(reactome, function(g) {
  length(intersect(graphite::nodes(g), row.names(CESC_EXP)))
  })
reactHuge_CESC <- reactome[nodesLength >= 10]

#Creating the multiomics
multiOmics_CESC_pathway <- Omics(data = list(expr = CESC_EXP, met = as.matrix(CESC_Methylation), mut = CESC_Mutation, cnv = CESC_CNV), methods = c("summarizeWithPca", "summarizeInCluster", "summarizeToBinaryEvents", "summarizeToBinaryEvents"), specificArgs = list(pcaArgs = list(name = "exp", shrink = "TRUE", method = "topological", maxPCs = 3), clusterArgs = list(name = "met"), binaryArgs = list(name = "mut", binaryClassMin = 10), binaryArgs = list(name = "cnv", binaryClassMin = 10)))

#Concidering just the genes which are present in expression data frame
genesToConsider_CESC <- row.names(CESC_EXP)

#Survival pathway test
startTime <- Sys.time()
multiOmicsFullHuge_CESC <- lapply(reactHuge_CESC, function(g) {
        print(g@title)
        set.seed(1234)
        fcl = multiOmicsSurvivalPathwayTest(multiOmics_CESC_pathway, g, survival_CESC, useThisGenes =genesToConsider_CESC)
        fcl
    })
endTime <- Sys.time()
print(endTime - startTime)
save(multiOmicsFullHuge_CESC, file = "./multiOmicsFullHuge_CESC.RData")

#Testing the result
multiOmicsFullHuge_CESC <- multiOmicsFullHuge_CESC[!duplicated(multiOmicsFullHuge_CESC)]
multiPathwayReport(multiOmicsFullHuge_CESC)

#Resampling step
pathwaySummaryHuge <- multiPathwayReport(multiOmicsFullHuge_CESC)
useThisPathwaysFM <- row.names(pathwaySummaryHuge[pathwaySummaryHuge$pvalue <= 0.05, ])
sPathwayM <- pathwaySummaryHuge[useThisPathwaysFM, , drop = T]
startTime <- Sys.time()
permsPathwaysHuge_CESC <- MOSClip:::resamplingPathway(fullMultiOmics = multiOmics_CESC_pathway, survival_CESC, reactHuge_CESC, nperm = 10, pathwaySubset = useThisPathwaysFM)
endTime <- Sys.time()
print(endTime - startTime)
save(permsPathwaysHuge_CESC, file = "./permsPathwaysHuge_CESC.RData")

#Adding the resampling score to the previous matrix
pathwaySummaryHuge <- multiPathwayReport(multiOmicsFullHuge_CESC)
useThisPathwaysFM <- row.names(pathwaySummaryHuge[pathwaySummaryHuge$pvalue <= 0.05, ])
sPathwayM <- pathwaySummaryHuge[useThisPathwaysFM, , drop = T]
resampligSuccessCount <- getPathwaysModulesSuccess(perms = permsPathwaysHuge_CESC, moduleSummary = sPathwayM)
pathwaySummaryHuge_CESC <- addResamplingCounts(pathwaySummaryHuge, resampligSuccessCount)

#Chossing those pathways with a score of at least 8
pathwaySummaryHuge_CESC <- pathwaySummaryHuge_CESC[pathwaySummaryHuge_CESC$resamplingCount >= 8,]
save(pathwaySummaryHuge_CESC, file = "./pathwaySummaryHuge_CESC.RData")
pathwaySummaryHuge_CESC

### 2.3 OV 

#Filtering the reactome pathways
#choosing only the ones who  have  expression
nodesLength <- sapply(reactome, function(g) {
  length(intersect(graphite::nodes(g), row.names(OV_EXP)))
  })
reactHuge_OV <- reactome[nodesLength >= 10]

#Creating the multiomics
multiOmics_OV_pathway <- Omics(data = list(expr = OV_EXP, met = as.matrix(OV_Methylation), mut = OV_Mutation, cnv = OV_CNV), methods = c("summarizeWithPca", "summarizeInCluster", "summarizeToBinaryEvents", "summarizeToBinaryEvents"), specificArgs = list(pcaArgs = list(name = "exp", shrink = "TRUE", method = "topological", maxPCs = 3), clusterArgs = list(name = "met"), binaryArgs = list(name = "mut", binaryClassMin = 10), binaryArgs = list(name = "cnv", binaryClassMin = 10)))

#Concidering just the genes which are present in expression data frame
genesToConsider_OV <- row.names(OV_EXP)

#Survival pathway test
startTime <- Sys.time()
multiOmicsFullHuge_OV <- lapply(reactHuge_OV, function(g) {
        print(g@title)
        set.seed(1234)
        fcl = multiOmicsSurvivalPathwayTest(multiOmics_OV_pathway, g, survival_OV, useThisGenes =genesToConsider_OV)
        fcl
    })
endTime <- Sys.time()
print(endTime - startTime)
save(multiOmicsFullHuge_OV, file = "./multiOmicsFullHuge_OV.RData")

#Testing the result
multiOmicsFullHuge_OV <- multiOmicsFullHuge_OV[!duplicated(multiOmicsFullHuge_OV)]
multiPathwayReport(multiOmicsFullHuge_OV)

#Resampling step
pathwaySummaryHuge <- multiPathwayReport(multiOmicsFullHuge_OV)
useThisPathwaysFM <- row.names(pathwaySummaryHuge[pathwaySummaryHuge$pvalue <= 0.05, ])
sPathwayM <- pathwaySummaryHuge[useThisPathwaysFM, , drop = T]
startTime <- Sys.time()
permsPathwaysHuge_OV <- MOSClip:::resamplingPathway(fullMultiOmics = multiOmics_OV_pathway, survival_OV, reactHuge_OV, nperm = 10, pathwaySubset = useThisPathwaysFM)
endTime <- Sys.time()
print(endTime - startTime)
save(permsPathwaysHuge_OV, file = "./permsPathwaysHuge_OV.RData")

#Adding the resampling score to the previous matrix
pathwaySummaryHuge <- multiPathwayReport(multiOmicsFullHuge_OV)
useThisPathwaysFM <- row.names(pathwaySummaryHuge[pathwaySummaryHuge$pvalue <= 0.05, ])
sPathwayM <- pathwaySummaryHuge[useThisPathwaysFM, , drop = T]
resampligSuccessCount <- getPathwaysModulesSuccess(perms = permsPathwaysHuge_OV, moduleSummary = sPathwayM)
pathwaySummaryHuge_OV <- addResamplingCounts(pathwaySummaryHuge, resampligSuccessCount)

#Chossing those pathways with a score of at least 8
pathwaySummaryHuge_OV <- pathwaySummaryHuge_OV[pathwaySummaryHuge_OV$resamplingCount >= 8,]
save(pathwaySummaryHuge_OV, file = "./pathwaySummaryHuge_OV.RData")
pathwaySummaryHuge_OV

### 2.4 UCEC

#Filtering the reactome pathways
#choosing only the ones who  have  expression
nodesLength <- sapply(reactome, function(g) {
  length(intersect(graphite::nodes(g), row.names(UCEC_EXP)))
  })
reactHuge_UCEC <- reactome[nodesLength >= 10]

#Creating the multiomics
multiOmics_UCEC_pathway <- Omics(data = list(expr = UCEC_EXP, met = as.matrix(UCEC_Methylation), mut = UCEC_Mutation, cnv = UCEC_CNV), methods = c("summarizeWithPca", "summarizeInCluster", "summarizeToBinaryEvents", "summarizeToBinaryEvents"), specificArgs = list(pcaArgs = list(name = "exp", shrink = "TRUE", method = "topological", maxPCs = 3), clusterArgs = list(name = "met"), binaryArgs = list(name = "mut", binaryClassMin = 10), binaryArgs = list(name = "cnv", binaryClassMin = 10)))

#Concidering just the genes which are present in expression data frame
genesToConsider_UCEC <- row.names(UCEC_EXP)

#Survival pathway test
startTime <- Sys.time()
multiOmicsFullHuge_UCEC <- lapply(reactHuge_UCEC, function(g) {
        print(g@title)
        set.seed(1234)
        fcl = multiOmicsSurvivalPathwayTest(multiOmics_UCEC_pathway, g, survival_UCEC, useThisGenes =genesToConsider_UCEC)
        fcl
    })
endTime <- Sys.time()
print(endTime - startTime)
save(multiOmicsFullHuge_UCEC, file = "./multiOmicsFullHuge_UCEC.RData")

#Testing the result
multiOmicsFullHuge_UCEC <- multiOmicsFullHuge_UCEC[!duplicated(multiOmicsFullHuge_UCEC)]
multiPathwayReport(multiOmicsFullHuge_UCEC)

#Resampling step
pathwaySummaryHuge <- multiPathwayReport(multiOmicsFullHuge_UCEC)
useThisPathwaysFM <- row.names(pathwaySummaryHuge[pathwaySummaryHuge$pvalue <= 0.05, ])
sPathwayM <- pathwaySummaryHuge[useThisPathwaysFM, , drop = T]
startTime <- Sys.time()
permsPathwaysHuge_UCEC <- MOSClip:::resamplingPathway(fullMultiOmics = multiOmics_UCEC_pathway, survival_UCEC, reactHuge_UCEC, nperm = 10, pathwaySubset = useThisPathwaysFM)
endTime <- Sys.time()
print(endTime - startTime)
save(permsPathwaysHuge_UCEC, file = "./permsPathwaysHuge_UCEC.RData")

#Adding the resampling score to the previous matrix
pathwaySummaryHuge <- multiPathwayReport(multiOmicsFullHuge_UCEC)
useThisPathwaysFM <- row.names(pathwaySummaryHuge[pathwaySummaryHuge$pvalue <= 0.05, ])
sPathwayM <- pathwaySummaryHuge[useThisPathwaysFM, , drop = T]
resampligSuccessCount <- getPathwaysModulesSuccess(perms = permsPathwaysHuge_UCEC, moduleSummary = sPathwayM)
pathwaySummaryHuge_UCEC <- addResamplingCounts(pathwaySummaryHuge, resampligSuccessCount)

#Chossing those pathways with a score of at least 8
pathwaySummaryHuge_UCEC <- pathwaySummaryHuge_UCEC[pathwaySummaryHuge_UCEC$resamplingCount >= 8,]
save(pathwaySummaryHuge_UCEC, file = "./pathwaySummaryHuge_UCEC.RData")
pathwaySummaryHuge_UCEC

### 2.5 UCS 

#Filtering the reactome pathways
#choosing only the ones who  have  expression
nodesLength <- sapply(reactome, function(g) {
  length(intersect(graphite::nodes(g), row.names(UCS_EXP)))
  })
reactHuge_UCS <- reactome[nodesLength >= 10]

#Creating the multiomics
multiOmics_UCS_pathway <- Omics(data = list(expr = UCS_EXP, met = as.matrix(UCS_Methylation), mut = UCS_Mutation, cnv = UCS_CNV), methods = c("summarizeWithPca", "summarizeInCluster", "summarizeToBinaryEvents", "summarizeToBinaryEvents"), specificArgs = list(pcaArgs = list(name = "exp", shrink = "TRUE", method = "topological", maxPCs = 3), clusterArgs = list(name = "met"), binaryArgs = list(name = "mut", binaryClassMin = 10), binaryArgs = list(name = "cnv", binaryClassMin = 10)))

#Concidering just the genes which are present in expression data frame
genesToConsider_UCS <- row.names(UCS_EXP)

#Survival pathway test
startTime <- Sys.time()
multiOmicsFullHuge_UCS <- lapply(reactHuge_UCS, function(g) {
        print(g@title)
        set.seed(1234)
        fcl = multiOmicsSurvivalPathwayTest(multiOmics_UCS_pathway, g, survival_UCS, useThisGenes =genesToConsider_UCS)
        fcl
    })
endTime <- Sys.time()
print(endTime - startTime)
save(multiOmicsFullHuge_UCS, file = "./multiOmicsFullHuge_UCS.RData")

#Testing the result
multiOmicsFullHuge_UCS <- multiOmicsFullHuge_UCS[!duplicated(multiOmicsFullHuge_UCS)]
multiPathwayReport(multiOmicsFullHuge_UCS)

#Resampling step
pathwaySummaryHuge <- multiPathwayReport(multiOmicsFullHuge_UCS)
useThisPathwaysFM <- row.names(pathwaySummaryHuge[pathwaySummaryHuge$pvalue <= 0.05, ])
sPathwayM <- pathwaySummaryHuge[useThisPathwaysFM, , drop = T]
startTime <- Sys.time()
permsPathwaysHuge_UCS <- MOSClip:::resamplingPathway(fullMultiOmics = multiOmics_UCS_pathway, survival_UCS, reactHuge_UCS, nperm = 10, pathwaySubset = useThisPathwaysFM)
endTime <- Sys.time()
print(endTime - startTime)
save(permsPathwaysHuge_UCS, file = "./permsPathwaysHuge_UCS.RData")

#Adding the resampling score to the previous matrix
pathwaySummaryHuge <- multiPathwayReport(multiOmicsFullHuge_UCS)
useThisPathwaysFM <- row.names(pathwaySummaryHuge[pathwaySummaryHuge$pvalue <= 0.05, ])
sPathwayM <- pathwaySummaryHuge[useThisPathwaysFM, , drop = T]
resampligSuccessCount <- getPathwaysModulesSuccess(perms = permsPathwaysHuge_UCS, moduleSummary = sPathwayM)
pathwaySummaryHuge_UCS <- addResamplingCounts(pathwaySummaryHuge, resampligSuccessCount)

#Chossing those pathways with a score of at least 8
pathwaySummaryHuge_UCS <- pathwaySummaryHuge_UCS[pathwaySummaryHuge_UCS$resamplingCount >= 8,]
save(pathwaySummaryHuge_UCS, file = "./pathwaySummaryHuge_UCS.RData")
pathwaySummaryHuge_UCS

## 3. MOSClip module test

### 3.1 BRCA

#Making the reactHuge:
nodesLength <- sapply(reactome, function(g) {
  length(intersect(graphite::nodes(g), row.names(BRCA_EXP)))
  })
reactHuge_BRCA <- reactome[nodesLength >= 10]

#Performing the module analysis only on the significant pathways
pathway_list_BRCA <- list()
for (i in rownames(pathwaySummaryHuge_BRCA)){
  pathway_list_BRCA <- c(pathway_list_BRCA, reactHuge_BRCA@entries[i])
}

genesToConsider <- row.names(BRCA_EXP)
multiOmics_BRCA_modules <- Omics(data = list(expr = BRCA_EXP, met = as.matrix(BRCA_Methylation), mut = BRCA_Mutation, cnv = BRCA_CNV), methods = c("summarizeWithPca", "summarizeInCluster", "summarizeToBinaryEvents", "summarizeToBinaryEvents"), specificArgs = list(pcaArgs = list(name = "exp", shrink = "FALSE", method = "sparse", maxPCs = 3), clusterArgs = list(name = "met"), binaryArgs = list(name = "mut", binaryClassMin = 10), binaryArgs = list(name = "cnv", binaryClassMin = 10)))
startTime <- Sys.time()
multiOmics_module_BRCA <- lapply(pathway_list_BRCA, function(g) {
        print(g@title) 
        set.seed(1234)
        fcl = multiOmicsSurvivalModuleTest(multiOmics_BRCA_modules, g, survival_BRCA, useThisGenes = genesToConsider)
        fcl
    })
endTime <- Sys.time()
print(endTime - startTime)
save(multiOmics_module_BRCA, file = "./multiOmics_module_BRCA.RData")

#resampling
moduleSummary <- multiPathwayModuleReport(multiOmics_module_BRCA)
useThisPathways <- unique(moduleSummary$pathway[moduleSummary$pvalue <= 0.05])
sModule <- moduleSummary[moduleSummary$pathway %in% useThisPathways, , drop = T]
startTime <- Sys.time()
permsModules_BRCA <- resampling(fullMultiOmics = multiOmics_BRCA_modules, survival_BRCA, pathway_list_BRCA, nperm = 10, pathwaySubset = useThisPathways)
endTime <- Sys.time()
print(endTime - startTime)
save(permsModules_BRCA, file = "./permsModules_BRCA.RData")

#module analysis
moduleSummary <- multiPathwayModuleReport(multiOmics_module_BRCA)
useThisPathways <- unique(moduleSummary$pathway[moduleSummary$pvalue <= 0.05])
sModule <- moduleSummary[moduleSummary$pathway %in% useThisPathways, , drop = T]
resampligSuccessCount_module <- getPathwaysModulesSuccess(perms = permsModules_BRCA, moduleSummary = sModule)
moduleSummary <- multiPathwayModuleReport(multiOmics_module_BRCA)
moduleSummary_BRCA <- addResamplingCounts(moduleSummary, resampligSuccessCount_module)
moduleSummary_BRCA <- moduleSummary_BRCA[!is.na(moduleSummary_BRCA$pathway),]

#Choosing modules with p-value of less than or equal to 0.05 and resampling score of at least 8
moduleSummary_BRCA <- moduleSummary_BRCA[moduleSummary_BRCA$pvalue <= 0.05 & moduleSummary_BRCA$resamplingCount >= 8,]

#Correcting the p-value
p_values <- p.adjust(moduleSummary_BRCA$pvalue)
moduleSummary_BRCA$pvalue <- p_values
moduleSummary_BRCA <- moduleSummary_BRCA[moduleSummary_BRCA$pvalue <= 0.1, ]
moduleSummary_BRCA
save(moduleSummary_BRCA, file = "./ModuleSummary_BRCA.RData")

### 3.2 CESC

#Making the reactHuge
nodesLength <- sapply(reactome, function(g) {
  length(intersect(graphite::nodes(g), row.names(CESC_EXP)))
  })
reactHuge_CESC <- reactome[nodesLength >= 10]

#Performing the module analysis only on the significant pathways
pathway_list_CESC <- list()
for (i in rownames(pathwaySummaryHuge_CESC)){
  pathway_list_CESC <- c(pathway_list_CESC, reactHuge_CESC@entries[i])
}

genesToConsider <- row.names(CESC_EXP)
multiOmics_CESC_modules <- Omics(data = list(expr = CESC_EXP, met = as.matrix(CESC_Methylation), mut = CESC_Mutation, cnv = CESC_CNV), methods = c("summarizeWithPca", "summarizeInCluster", "summarizeToBinaryEvents", "summarizeToBinaryEvents"), specificArgs = list(pcaArgs = list(name = "exp", shrink = "FALSE", method = "sparse", maxPCs = 3), clusterArgs = list(name = "met"), binaryArgs = list(name = "mut", binaryClassMin = 10), binaryArgs = list(name = "cnv", binaryClassMin = 10)))
startTime <- Sys.time()
multiOmics_module_CESC <- lapply(pathway_list_CESC, function(g) {
        print(g@title) 
        set.seed(1234)
        fcl = multiOmicsSurvivalModuleTest(multiOmics_CESC_modules, g, survival_CESC, useThisGenes = genesToConsider)
        fcl
    })
endTime <- Sys.time()
print(endTime - startTime)
save(multiOmics_module_CESC, file = "./multiOmics_module_CESC.RData")

#resampling
moduleSummary <- multiPathwayModuleReport(multiOmics_module_CESC)
useThisPathways <- unique(moduleSummary$pathway[moduleSummary$pvalue <= 0.05])
sModule <- moduleSummary[moduleSummary$pathway %in% useThisPathways, , drop = T]
startTime <- Sys.time()
permsModules_CESC <- resampling(fullMultiOmics = multiOmics_CESC_modules, survival_CESC, pathway_list_CESC, nperm = 10, pathwaySubset = useThisPathways)
endTime <- Sys.time()
print(endTime - startTime)
save(permsModules_CESC, file = "./permsModules_CESC.RData")

#module analysis
moduleSummary <- multiPathwayModuleReport(multiOmics_module_CESC)
useThisPathways <- unique(moduleSummary$pathway[moduleSummary$pvalue <= 0.05])
sModule <- moduleSummary[moduleSummary$pathway %in% useThisPathways, , drop = T]
resampligSuccessCount_module <- getPathwaysModulesSuccess(perms = permsModules_CESC, moduleSummary = sModule)
moduleSummary <- multiPathwayModuleReport(multiOmics_module_CESC)
moduleSummary_CESC <- addResamplingCounts(moduleSummary, resampligSuccessCount_module)
moduleSummary_CESC <- moduleSummary_CESC[!is.na(moduleSummary_CESC$pathway),]

#Choosing modules with p-value of less than or equal to 0.05 and resampling score of at least 8
moduleSummary_CESC <- moduleSummary_CESC[moduleSummary_CESC$pvalue <= 0.05 & moduleSummary_CESC$resamplingCount >= 8,]

#Correcting the p-value
p_values <- p.adjust(moduleSummary_CESC$pvalue)
moduleSummary_CESC$pvalue <- p_values
moduleSummary_CESC <- moduleSummary_CESC[moduleSummary_CESC$pvalue <= 0.1, ]
moduleSummary_CESC
save(moduleSummary_CESC, file = "./ModuleSummary_CESC.RData")

### 3.3 OV

#Making the reactHuge
nodesLength <- sapply(reactome, function(g) {
  length(intersect(graphite::nodes(g), row.names(OV_EXP)))
  })
reactHuge_OV <- reactome[nodesLength >= 10]

#Performing the module analysis only on the significant pathways
pathway_list_OV <- list()
for (i in rownames(pathwaySummaryHuge_OV)){
  pathway_list_OV <- c(pathway_list_OV, reactHuge_OV@entries[i])
}

genesToConsider <- row.names(OV_EXP)
multiOmics_OV_modules <- Omics(data = list(expr = OV_EXP, met = as.matrix(OV_Methylation), mut = OV_Mutation, cnv = OV_CNV), methods = c("summarizeWithPca", "summarizeInCluster", "summarizeToBinaryEvents", "summarizeToBinaryEvents"), specificArgs = list(pcaArgs = list(name = "exp", shrink = "FALSE", method = "sparse", maxPCs = 3), clusterArgs = list(name = "met"), binaryArgs = list(name = "mut", binaryClassMin = 10), binaryArgs = list(name = "cnv", binaryClassMin = 10)))
startTime <- Sys.time()
multiOmics_module_OV <- lapply(pathway_list_OV, function(g) {
        print(g@title) 
        set.seed(1234)
        fcl = multiOmicsSurvivalModuleTest(multiOmics_OV_modules, g, survival_OV, useThisGenes = genesToConsider)
        fcl
    })
endTime <- Sys.time()
print(endTime - startTime)
save(multiOmics_module_OV, file = "./multiOmics_module_OV.RData")

#resampling
moduleSummary <- multiPathwayModuleReport(multiOmics_module_OV)
useThisPathways <- unique(moduleSummary$pathway[moduleSummary$pvalue <= 0.05])
sModule <- moduleSummary[moduleSummary$pathway %in% useThisPathways, , drop = T]
startTime <- Sys.time()
permsModules_OV <- resampling(fullMultiOmics = multiOmics_OV_modules, survival_OV, pathway_list_OV, nperm = 10, pathwaySubset = useThisPathways)
endTime <- Sys.time()
print(endTime - startTime)
save(permsModules_OV, file = "./permsModules_OV.RData")

#module analysis
moduleSummary <- multiPathwayModuleReport(multiOmics_module_OV)
useThisPathways <- unique(moduleSummary$pathway[moduleSummary$pvalue <= 0.05])
sModule <- moduleSummary[moduleSummary$pathway %in% useThisPathways, , drop = T]
resampligSuccessCount_module <- getPathwaysModulesSuccess(perms = permsModules_OV, moduleSummary = sModule)
moduleSummary <- multiPathwayModuleReport(multiOmics_module_OV)
moduleSummary_OV <- addResamplingCounts(moduleSummary, resampligSuccessCount_module)
moduleSummary_OV <- moduleSummary_OV[!is.na(moduleSummary_OV$pathway),]

#Choosing modules with p-value of less than or equal to 0.05 and resampling score of at least 8
moduleSummary_OV <- moduleSummary_OV[moduleSummary_OV$pvalue <= 0.05 & moduleSummary_OV$resamplingCount >= 8,]

#Correcting the p-value
p_values <- p.adjust(moduleSummary_OV$pvalue)
moduleSummary_OV$pvalue <- p_values
moduleSummary_OV <- moduleSummary_OV[moduleSummary_OV$pvalue <= 0.1, ]
moduleSummary_OV
save(moduleSummary_OV, file = "./ModuleSummary_OV.RData")

### 3.4 UCEC

#Making the reactHuge
nodesLength <- sapply(reactome, function(g) {
  length(intersect(graphite::nodes(g), row.names(UCEC_EXP)))
  })
reactHuge_UCEC <- reactome[nodesLength >= 10]

#Performing the module analysis only on the significant pathways
pathway_list_UCEC <- list()
for (i in rownames(pathwaySummaryHuge_UCEC)){
  pathway_list_UCEC <- c(pathway_list_UCEC, reactHuge_UCEC@entries[i])
}

genesToConsider <- row.names(UCEC_EXP)
multiOmics_UCEC_modules <- Omics(data = list(expr = UCEC_EXP, met = as.matrix(UCEC_Methylation), mut = UCEC_Mutation, cnv = UCEC_CNV), methods = c("summarizeWithPca", "summarizeInCluster", "summarizeToBinaryEvents", "summarizeToBinaryEvents"), specificArgs = list(pcaArgs = list(name = "exp", shrink = "FALSE", method = "sparse", maxPCs = 3), clusterArgs = list(name = "met"), binaryArgs = list(name = "mut", binaryClassMin = 10), binaryArgs = list(name = "cnv", binaryClassMin = 10)))
startTime <- Sys.time()
multiOmics_module_UCEC <- lapply(pathway_list_UCEC, function(g) {
        print(g@title) 
        set.seed(1234)
        fcl = multiOmicsSurvivalModuleTest(multiOmics_UCEC_modules, g, survival_UCEC, useThisGenes = genesToConsider)
        fcl
    })
endTime <- Sys.time()
print(endTime - startTime)
save(multiOmics_module_UCEC, file = "./multiOmics_module_UCEC.RData")

#resampling
moduleSummary <- multiPathwayModuleReport(multiOmics_module_UCEC)
useThisPathways <- unique(moduleSummary$pathway[moduleSummary$pvalue <= 0.05])
sModule <- moduleSummary[moduleSummary$pathway %in% useThisPathways, , drop = T]
startTime <- Sys.time()
permsModules_UCEC <- resampling(fullMultiOmics = multiOmics_UCEC_modules, survival_UCEC, pathway_list_UCEC, nperm = 10, pathwaySubset = useThisPathways)
endTime <- Sys.time()
print(endTime - startTime)
save(permsModules_UCEC, file = "./permsModules_UCEC.RData")

#module analysis
moduleSummary <- multiPathwayModuleReport(multiOmics_module_UCEC)
useThisPathways <- unique(moduleSummary$pathway[moduleSummary$pvalue <= 0.05])
sModule <- moduleSummary[moduleSummary$pathway %in% useThisPathways, , drop = T]
resampligSuccessCount_module <- getPathwaysModulesSuccess(perms = permsModules_UCEC, moduleSummary = sModule)
moduleSummary <- multiPathwayModuleReport(multiOmics_module_UCEC)
moduleSummary_UCEC <- addResamplingCounts(moduleSummary, resampligSuccessCount_module)
moduleSummary_UCEC <- moduleSummary_UCEC[!is.na(moduleSummary_UCEC$pathway),]

#Choosing modules with p-value of less than or equal to 0.05 and resampling score of at least 8
moduleSummary_UCEC <- moduleSummary_UCEC[moduleSummary_UCEC$pvalue <= 0.05 & moduleSummary_UCEC$resamplingCount >= 8,]

#Correcting the p-value
p_values <- p.adjust(moduleSummary_UCEC$pvalue)
moduleSummary_UCEC$pvalue <- p_values
moduleSummary_UCEC <- moduleSummary_UCEC[moduleSummary_UCEC$pvalue <= 0.1, ]
moduleSummary_UCEC
save(moduleSummary_UCEC, file = "./ModuleSummary_UCEC.RData")

### 3.5 UCS

#Making the reactHuge
nodesLength <- sapply(reactome, function(g) {
  length(intersect(graphite::nodes(g), row.names(UCS_EXP)))
  })
reactHuge_UCS <- reactome[nodesLength >= 10]

#Performing the module analysis only on the significant pathways
pathway_list_UCS <- list()
for (i in rownames(pathwaySummaryHuge_UCS)){
  pathway_list_UCS <- c(pathway_list_UCS, reactHuge_UCS@entries[i])
}

genesToConsider <- row.names(UCS_EXP)
multiOmics_UCS_modules <- Omics(data = list(expr = UCS_EXP, met = as.matrix(UCS_Methylation), mut = UCS_Mutation, cnv = UCS_CNV), methods = c("summarizeWithPca", "summarizeInCluster", "summarizeToBinaryEvents", "summarizeToBinaryEvents"), specificArgs = list(pcaArgs = list(name = "exp", shrink = "FALSE", method = "sparse", maxPCs = 3), clusterArgs = list(name = "met"), binaryArgs = list(name = "mut", binaryClassMin = 10), binaryArgs = list(name = "cnv", binaryClassMin = 10)))
startTime <- Sys.time()
multiOmics_module_UCS <- lapply(pathway_list_UCS, function(g) {
        print(g@title) 
        set.seed(1234)
        fcl = multiOmicsSurvivalModuleTest(multiOmics_UCS_modules, g, survival_UCS, useThisGenes = genesToConsider)
        fcl
    })
endTime <- Sys.time()
print(endTime - startTime)
save(multiOmics_module_UCS, file = "./multiOmics_module_UCS.RData")

#resampling
moduleSummary <- multiPathwayModuleReport(multiOmics_module_UCS)
useThisPathways <- unique(moduleSummary$pathway[moduleSummary$pvalue <= 0.05])
sModule <- moduleSummary[moduleSummary$pathway %in% useThisPathways, , drop = T]
startTime <- Sys.time()
permsModules_UCS <- resampling(fullMultiOmics = multiOmics_UCS_modules, survival_UCS, pathway_list_UCS, nperm = 10, pathwaySubset = useThisPathways)
endTime <- Sys.time()
print(endTime - startTime)
save(permsModules_UCS, file = "./permsModules_UCS.RData")

#module analysis
moduleSummary <- multiPathwayModuleReport(multiOmics_module_UCS)
useThisPathways <- unique(moduleSummary$pathway[moduleSummary$pvalue <= 0.05])
sModule <- moduleSummary[moduleSummary$pathway %in% useThisPathways, , drop = T]
resampligSuccessCount_module <- getPathwaysModulesSuccess(perms = permsModules_UCS, moduleSummary = sModule)
moduleSummary <- multiPathwayModuleReport(multiOmics_module_UCS)
moduleSummary_UCS <- addResamplingCounts(moduleSummary, resampligSuccessCount_module)
moduleSummary_UCS <- moduleSummary_UCS[!is.na(moduleSummary_UCS$pathway),]

#Choosing modules with p-value of less than or equal to 0.05 and resampling score of at least 8
moduleSummary_UCS <- moduleSummary_UCS[moduleSummary_UCS$pvalue <= 0.05 & moduleSummary_UCS$resamplingCount >= 8,]

#Correcting the p-value
p_values <- p.adjust(moduleSummary_UCS$pvalue)
moduleSummary_UCS$pvalue <- p_values
moduleSummary_UCS <- moduleSummary_UCS[moduleSummary_UCS$pvalue <= 0.15, ]
moduleSummary_UCS
save(moduleSummary_UCS, file = "./ModuleSummary_UCS.RData")

